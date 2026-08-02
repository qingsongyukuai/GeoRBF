use crate::functional::CanonicalFunctional;
use crate::math::{canonical_zero, dot3};
use crate::numerical::EQUALITY_KKT_POLICY_V1;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GlobalAnisotropyMetric {
    matrix: [[f64; 3]; 3],
}

impl GlobalAnisotropyMetric {
    pub(crate) fn new(matrix: [[f64; 3]; 3]) -> Result<Self, MetricError> {
        for (row, values) in matrix.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                if !value.is_finite() {
                    return Err(MetricError::NonFinite { row, column });
                }
            }
        }
        for row in 0..3 {
            for column in (row + 1)..3 {
                if matrix[row][column] != matrix[column][row] {
                    return Err(MetricError::NotSymmetric { row, column });
                }
            }
        }

        let leading_two = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        let determinant = determinant(matrix);
        if !determinant.is_finite() {
            return Err(MetricError::NonFiniteDeterminant);
        }
        if matrix[0][0] <= 0.0 || leading_two <= 0.0 || determinant <= 0.0 {
            return Err(MetricError::NotPositiveDefinite);
        }
        let determinant_tolerance = EQUALITY_KKT_POLICY_V1.metric_determinant_one_multiplier
            * f64::EPSILON
            * determinant.abs().max(1.0);
        if (determinant - 1.0).abs() > determinant_tolerance {
            return Err(MetricError::DeterminantNotOne { determinant });
        }
        Ok(Self { matrix })
    }

    pub(crate) fn identity() -> Self {
        Self {
            matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    pub(crate) fn matrix(&self) -> [[f64; 3]; 3] {
        self.matrix
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum MetricError {
    NonFinite { row: usize, column: usize },
    NotSymmetric { row: usize, column: usize },
    NotPositiveDefinite,
    NonFiniteDeterminant,
    DeterminantNotOne { determinant: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CubicJet {
    value: f64,
    gradient_x: [f64; 3],
    gradient_y: [f64; 3],
    mixed_xy: [[f64; 3]; 3],
}

impl CubicJet {
    pub(crate) fn value(self) -> f64 {
        self.value
    }

    pub(crate) fn gradient_x(self) -> [f64; 3] {
        self.gradient_x
    }

    pub(crate) fn gradient_y(self) -> [f64; 3] {
        self.gradient_y
    }

    pub(crate) fn mixed_xy(self) -> [[f64; 3]; 3] {
        self.mixed_xy
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CubicKernel;

impl CubicKernel {
    pub(crate) fn jet(x: [f64; 3], y: [f64; 3], metric: &GlobalAnisotropyMetric) -> CubicJet {
        let delta = [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
        let displacement_scale = delta
            .iter()
            .map(|component| component.abs())
            .fold(0.0_f64, f64::max);
        if displacement_scale == 0.0 {
            return CubicJet {
                value: 0.0,
                gradient_x: [0.0; 3],
                gradient_y: [0.0; 3],
                mixed_xy: [[0.0; 3]; 3],
            };
        }
        let scaled_delta = delta.map(|component| component / displacement_scale);
        let scaled_metric_delta = multiply(metric.matrix, scaled_delta);
        let scaled_radius = dot3(scaled_delta, scaled_metric_delta).sqrt();
        let radius = displacement_scale * scaled_radius;
        let gradient_scale = 3.0 * displacement_scale * displacement_scale * scaled_radius;
        let gradient_x = scaled_metric_delta.map(|value| canonical_zero(gradient_scale * value));
        let gradient_y = gradient_x.map(|value| canonical_zero(-value));
        let mixed_xy = std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                canonical_zero(
                    -3.0 * displacement_scale
                        * (scaled_radius * metric.matrix[row][column]
                            + scaled_metric_delta[row] * scaled_metric_delta[column]
                                / scaled_radius),
                )
            })
        });
        CubicJet {
            value: radius * radius * radius,
            gradient_x,
            gradient_y,
            mixed_xy,
        }
    }

    pub(crate) fn pairing(
        left: &CanonicalFunctional,
        right: &CanonicalFunctional,
        metric: &GlobalAnisotropyMetric,
    ) -> f64 {
        left.terms()
            .iter()
            .flat_map(|left_term| {
                right.terms().iter().map(move |right_term| {
                    let jet = Self::jet(left_term.support(), right_term.support(), metric);
                    let left_gradient = left_term.gradient_coefficient();
                    let right_gradient = right_term.gradient_coefficient();
                    left_term.value_coefficient() * right_term.value_coefficient() * jet.value
                        + left_term.value_coefficient() * dot3(right_gradient, jet.gradient_y)
                        + right_term.value_coefficient() * dot3(left_gradient, jet.gradient_x)
                        + dot3(left_gradient, multiply(jet.mixed_xy, right_gradient))
                })
            })
            .sum()
    }
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}

fn multiply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::functional::{CanonicalFunctional, FunctionalDimension, FunctionalTerm};
    use crate::oracle_fixture::{
        hex_scalar, hex_values, verify_artifact_identity, verify_fixture_identity,
    };

    const GENERAL_JET: &str =
        include_str!("../validation/oracle/cubic-v1/fixtures/cubic-general-jet.json");
    const ORIGIN: &str = include_str!("../validation/oracle/cubic-v1/fixtures/cubic-origin.json");
    const GENERALIZED_FUNCTIONAL: &str =
        include_str!("../validation/oracle/cubic-v1/fixtures/cubic-generalized-functional.json");
    const MANIFEST: &str = include_str!("../validation/oracle/cubic-v1/source-manifest.json");
    const GENERAL_JET_CASE: &str =
        include_str!("../validation/oracle/cubic-v1/cases/cubic-general-jet.json");
    const ORIGIN_CASE: &str = include_str!("../validation/oracle/cubic-v1/cases/cubic-origin.json");
    const GENERALIZED_FUNCTIONAL_CASE: &str =
        include_str!("../validation/oracle/cubic-v1/cases/cubic-generalized-functional.json");

    fn assert_close(actual: f64, expected: f64, ulps: f64) {
        let limit = ulps * f64::EPSILON * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= limit,
            "actual={actual:e}, expected={expected:e}, limit={limit:e}"
        );
    }

    fn metric() -> GlobalAnisotropyMetric {
        GlobalAnisotropyMetric::new([[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]])
            .expect("the oracle metric is SPD with determinant one")
    }

    #[test]
    fn global_metric_rejects_even_small_asymmetry_instead_of_repairing_it() {
        let failure =
            GlobalAnisotropyMetric::new([[1.0, 1.0e-15, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
                .expect_err("the resolved metric contract requires a symmetric matrix");

        assert_eq!(failure, MetricError::NotSymmetric { row: 0, column: 1 });
    }

    #[test]
    fn global_metric_rejects_a_nonfinite_computed_determinant() {
        let failure = GlobalAnisotropyMetric::new([
            [1.0e308, 0.0, 0.0],
            [0.0, 1.0e308, 0.0],
            [0.0, 0.0, 1.0e308],
        ])
        .expect_err("overflow cannot masquerade as a determinant-one metric");

        assert_eq!(failure, MetricError::NonFiniteDeterminant);
    }

    #[test]
    fn cubic_distinguishes_true_origin_from_squared_distance_underflow() {
        let jet = CubicKernel::jet(
            [1.0e-200, 0.0, 0.0],
            [0.0; 3],
            &GlobalAnisotropyMetric::identity(),
        );

        assert!((jet.mixed_xy()[0][0] + 6.0e-200).abs() <= 1.0e-214);
        assert!((jet.mixed_xy()[1][1] + 3.0e-200).abs() <= 1.0e-214);
        assert!((jet.mixed_xy()[2][2] + 3.0e-200).abs() <= 1.0e-214);
    }

    #[test]
    fn cubic_general_jet_matches_the_independent_fixture() {
        verify_fixture_identity(
            GENERAL_JET,
            "dbbe6bfbfe8836a253c1beac55b2662498a4775e4ad3f0cc6aa0f5c347e7ab46",
            "case-v1-fb096a2d18463c80697b610cf201a01cb895e752bcde547f7fe94d64e448a647",
            "GeoRBF issue #6",
            "sha256:a070035b6bface57abd2945358a0a9fa892dbe03400733108024f72974108368",
        );
        assert!(MANIFEST.contains("\"version\": \"risk-spike-15-v1\""));
        assert!(MANIFEST.contains("\"version\": \"3.12.3\""));
        assert!(
            MANIFEST.contains(
                "sha256:dbbe6bfbfe8836a253c1beac55b2662498a4775e4ad3f0cc6aa0f5c347e7ab46"
            )
        );
        verify_artifact_identity(
            MANIFEST,
            "ba516a0374c799e7e7ed3225d84093a2db1e8298e7dedc0e74569dc325e3ce54",
        );
        verify_artifact_identity(
            GENERAL_JET_CASE,
            "b24f2609b51c4cb646cc903478cabfb23e6f74d3e766d68b978939303ac9be26",
        );

        let jet = CubicKernel::jet([1.25, -2.0, 0.5], [-0.75, 1.5, 2.0], &metric());

        assert_close(jet.value(), hex_scalar(GENERAL_JET, "value"), 64.0);
        for (actual, expected) in
            jet.gradient_x()
                .into_iter()
                .zip(hex_values(GENERAL_JET, "gradient_x", 3))
        {
            assert_close(actual, expected, 256.0);
        }
        for (actual, expected) in
            jet.gradient_y()
                .into_iter()
                .zip(hex_values(GENERAL_JET, "gradient_y", 3))
        {
            assert_close(actual, expected, 256.0);
        }
        for (actual, expected) in
            jet.mixed_xy()
                .into_iter()
                .flatten()
                .zip(hex_values(GENERAL_JET, "mixed_xy", 9))
        {
            assert_close(actual, expected, 1024.0);
        }
    }

    #[test]
    fn cubic_origin_uses_the_exact_positive_zero_analytic_limit() {
        verify_fixture_identity(
            ORIGIN,
            "a43c91b2165c760cec7abfb05604f0253abdce41febdfb13f8c4544fe5928e97",
            "case-v1-d3562269dbc1bbbe0c94ed504452a03288396b7cfcac56c1f777e76628462086",
            "GeoRBF issue #6",
            "sha256:c7a0f2fa158836555c45a25eda226254529403831a0c377a4ad6f7c2d16d5e80",
        );
        verify_artifact_identity(
            ORIGIN_CASE,
            "0eae84c29abb6a7c71b1d97733f8cf69181426a6f83507c2444268c6ce179f14",
        );
        assert!(
            hex_values(ORIGIN, "result", 24)
                .into_iter()
                .all(|value| value.to_bits() == 0)
        );

        let jet = CubicKernel::jet([1.25, -2.0, 0.5], [1.25, -2.0, 0.5], &metric());
        let values = std::iter::once(jet.value())
            .chain(jet.gradient_x())
            .chain(jet.gradient_y())
            .chain(jet.mixed_xy().into_iter().flatten());
        assert!(values.into_iter().all(|value| value.to_bits() == 0));
    }

    #[test]
    fn generalized_functional_pairing_and_affine_observation_match_the_fixture() {
        verify_fixture_identity(
            GENERALIZED_FUNCTIONAL,
            "8ff330dd24c5bf6efb076ff4636da058451057713567ce3a77b2b86fb94bf6b0",
            "case-v1-ec9c5da89d9fa94f037a61ac791c3508def3238f96501aa9edae7b5673c86950",
            "GeoRBF issues #3, #6, and #10",
            "sha256:856035654212317863d49567179afd61a3ce30f8b5e1d2a256e47a305f14bbf1",
        );
        verify_artifact_identity(
            GENERALIZED_FUNCTIONAL_CASE,
            "a62164c189d26acf7d8d51f286e545a7aaf1d3a7a6bdc89bea3f2b34296128ad",
        );
        let left = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![FunctionalTerm::new([0.0; 3], 2.0, [1.0, -0.5, 0.25])],
        )
        .expect("left fixture functional is valid");
        let right = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![FunctionalTerm::new(
                [1.0, -2.0, 0.5],
                -1.5,
                [0.75, 1.0, -0.25],
            )],
        )
        .expect("right fixture functional is valid");

        assert_close(
            CubicKernel::pairing(&left, &right, &GlobalAnisotropyMetric::identity()),
            hex_scalar(GENERALIZED_FUNCTIONAL, "cubic_pairing"),
            1024.0,
        );
        let affine = [1.125, -0.75, 2.0, 0.5];
        let affine_observations = hex_values(
            GENERALIZED_FUNCTIONAL,
            "manufactured_affine_observations",
            2,
        );
        assert_eq!(
            left.evaluate_affine(affine[0], [affine[1], affine[2], affine[3]]),
            affine_observations[0]
        );
        assert_eq!(
            right.evaluate_affine(affine[0], [affine[1], affine[2], affine[3]]),
            affine_observations[1]
        );
    }

    #[test]
    fn cubic_exchange_derivative_signs_and_frame_metric_covariance_hold() {
        let x = [1.25, -2.0, 0.5];
        let y = [-0.75, 1.5, 2.0];
        let metric = metric();
        let forward = CubicKernel::jet(x, y, &metric);
        let exchanged = CubicKernel::jet(y, x, &metric);
        assert_eq!(forward.value(), exchanged.value());
        for axis in 0..3 {
            assert_eq!(forward.gradient_x()[axis], exchanged.gradient_y()[axis]);
            assert_eq!(forward.gradient_y()[axis], exchanged.gradient_x()[axis]);
            for other in 0..3 {
                assert_eq!(
                    forward.mixed_xy()[axis][other],
                    exchanged.mixed_xy()[other][axis]
                );
            }
        }

        let scale = 2.0;
        let translation = [3.0, -4.0, 1.0];
        let transform = |point: [f64; 3]| {
            [
                scale * point[1] + translation[0],
                scale * point[0] + translation[1],
                -scale * point[2] + translation[2],
            ]
        };
        let transformed_metric =
            GlobalAnisotropyMetric::new([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.5]])
                .expect("orthogonal metric covariance preserves determinant one");
        let transformed = CubicKernel::jet(transform(x), transform(y), &transformed_metric);
        assert_close(transformed.value(), scale.powi(3) * forward.value(), 64.0);
        let rotate = |vector: [f64; 3]| [vector[1], vector[0], -vector[2]];
        for (actual, expected) in transformed
            .gradient_x()
            .into_iter()
            .zip(rotate(forward.gradient_x()).map(|value| scale.powi(2) * value))
        {
            assert_close(actual, expected, 256.0);
        }
        let rotation = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]];
        for row in 0..3 {
            for column in 0..3 {
                let expected = (0..3)
                    .flat_map(|left| {
                        (0..3).map(move |right| {
                            rotation[row][left]
                                * forward.mixed_xy()[left][right]
                                * rotation[column][right]
                        })
                    })
                    .sum::<f64>()
                    * scale;
                assert_close(transformed.mixed_xy()[row][column], expected, 1024.0);
            }
        }
    }
}
