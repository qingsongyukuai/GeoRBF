use georbf::kernel::{FieldEnergyNormalization, FieldEnergyNormalizationError};
use georbf::observation::{
    FieldValueObservation, QuadraticPenalty, QuadraticPenaltyError, StandardDeviation,
    StandardDeviationError,
};
use georbf::problem::{BuildError, BuilderConfigurationError};
use georbf::problem::{FitConfiguration, ThreadBudget};
use georbf::{Point3, ProblemBuilder, SourceId};
use std::num::NonZeroUsize;

mod common {
    use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};

    pub fn frame() -> InputCoordinateFrame {
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .expect("the test frame is valid")
    }

    pub fn field_unit() -> FieldUnitLabel {
        FieldUnitLabel::new("stratigraphic-unit")
    }
}

#[test]
fn soft_field_value_scalars_accept_only_finite_positive_values() {
    let penalty = QuadraticPenalty::try_new(2.5).expect("a positive penalty is valid");
    let standard_deviation =
        StandardDeviation::try_new(0.25).expect("a positive standard deviation is valid");
    let normalization = FieldEnergyNormalization::try_new(4.0)
        .expect("a positive FieldEnergy normalization is valid");

    assert_eq!(penalty.weight(), 2.5);
    assert_eq!(standard_deviation.value(), 0.25);
    assert_eq!(normalization.factor(), 4.0);

    for invalid in [0.0, -1.0] {
        assert_eq!(
            QuadraticPenalty::try_new(invalid),
            Err(QuadraticPenaltyError::NotPositive)
        );
        assert_eq!(
            StandardDeviation::try_new(invalid),
            Err(StandardDeviationError::NotPositive)
        );
        assert_eq!(
            FieldEnergyNormalization::try_new(invalid),
            Err(FieldEnergyNormalizationError::NotPositive)
        );
    }

    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            QuadraticPenalty::try_new(invalid),
            Err(QuadraticPenaltyError::NotFinite)
        );
        assert_eq!(
            StandardDeviation::try_new(invalid),
            Err(StandardDeviationError::NotFinite)
        );
        assert_eq!(
            FieldEnergyNormalization::try_new(invalid),
            Err(FieldEnergyNormalizationError::NotFinite)
        );
    }
}

#[test]
fn soft_field_values_require_an_explicit_repairable_field_energy_normalization() {
    let mut problem = ProblemBuilder::new(common::frame(), common::field_unit());
    problem
        .add(
            FieldValueObservation::try_with_quadratic_penalty(
                SourceId::new("penalty"),
                Point3::try_new(0.0, 0.0, 0.0).unwrap(),
                1.5,
                QuadraticPenalty::try_new(2.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(
            FieldValueObservation::try_with_standard_deviation(
                SourceId::new("statistical"),
                Point3::try_new(1.0, 0.0, 0.0).unwrap(),
                2.0,
                StandardDeviation::try_new(0.5).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    problem.set_fit_configuration(
        FitConfiguration::default()
            .with_thread_budget(ThreadBudget::Exact(NonZeroUsize::new(2).unwrap())),
    );

    let failure = problem
        .build()
        .expect_err("a soft problem cannot hide its FieldEnergy scale");
    assert_eq!(
        failure.errors(),
        &[
            BuildError::MissingFieldEnergyNormalization,
            BuildError::UnsupportedThreadBudget { requested: 2 },
        ]
    );

    let mut problem = failure.into_builder();
    problem.set_fit_configuration(FitConfiguration::default());
    let normalization = FieldEnergyNormalization::try_new(3.0).unwrap();
    problem
        .set_field_energy_normalization(normalization)
        .expect("the retained builder can be repaired");
    assert_eq!(
        problem.set_field_energy_normalization(normalization),
        Err(BuilderConfigurationError::FieldEnergyNormalizationAlreadySet)
    );
    let snapshot = problem.build().expect("the repaired soft problem builds");

    assert_eq!(snapshot.observation_count(), 2);
    assert_eq!(snapshot.field_energy_normalization().factor(), 3.0);
}

fn manufactured_problem() -> ProblemBuilder {
    manufactured_problem_with_scales(1.0, 1.0)
}

fn manufactured_problem_with_scales(length_scale: f64, field_scale: f64) -> ProblemBuilder {
    let mut problem = ProblemBuilder::new(common::frame(), common::field_unit());
    let hard_values = [
        ([-1.0, -1.0, -1.0], 0.0),
        ([1.0, -1.0, -1.0], 1.0),
        ([-1.0, 1.0, -1.0], -0.5),
        ([-1.0, -1.0, 1.0], 1.5),
        ([1.0, 1.0, 0.5], 1.625),
    ];
    for (index, (support, value)) in hard_values.into_iter().enumerate() {
        let [x, y, z] = support.map(|coordinate| coordinate * length_scale);
        problem
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("hard-{index}")),
                    Point3::try_new(x, y, z).unwrap(),
                    field_scale * value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    problem
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected:.16e}, found {actual:.16e}, tolerance {tolerance:.3e}"
    );
}

#[test]
fn user_can_fit_a_soft_field_value_and_audit_its_physical_objective() {
    let soft_location = Point3::try_new(0.2, -0.3, 0.4).unwrap();
    let mut problem = manufactured_problem();
    problem
        .add(
            FieldValueObservation::try_with_quadratic_penalty(
                SourceId::new("soft-value"),
                soft_location,
                2.0,
                QuadraticPenalty::try_new(2.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();

    let fit = problem
        .build()
        .unwrap()
        .fit()
        .expect("the soft quadratic Equality problem should fit through faer");
    let report = fit.report();
    let assessment = &report.soft_field_values()[0];
    let sample = fit.model().evaluate(soft_location).unwrap();

    assert_eq!(report.problem_size().scalar_soft_relations(), 1);
    assert_eq!(report.problem_size().quadratic_objective_terms(), 1);
    assert_eq!(assessment.source_id(), &SourceId::new("soft-value"));
    assert_eq!(assessment.target(), 2.0);
    assert_eq!(assessment.quadratic_penalty().unwrap().weight(), 2.0);
    assert_eq!(assessment.standard_deviation(), None);
    assert_close(assessment.recovered_value(), sample.value(), 1.0e-11);
    assert_close(
        assessment.residual(),
        assessment.recovered_value() - assessment.target(),
        1.0e-14,
    );
    assert_close(
        assessment.loss(),
        assessment.residual() * assessment.residual(),
        1.0e-13,
    );
    assert_close(
        report.total_objective().unwrap(),
        0.5 * report.field_energy().unwrap() + assessment.loss(),
        1.0e-11,
    );
    assert!(report.canonical_acceptance().unwrap().objective_verified());

    let batch = fit
        .model()
        .evaluate_batch(&[soft_location, Point3::try_new(0.0, 0.0, 0.0).unwrap()])
        .unwrap();
    assert_eq!(batch[0], sample);
    assert_eq!(batch.len(), 2);
}

#[test]
fn standard_deviation_matches_equivalent_quadratic_precision_without_losing_semantics() {
    let soft_location = Point3::try_new(0.2, -0.3, 0.4).unwrap();

    let mut penalty_problem = manufactured_problem();
    penalty_problem
        .add(
            FieldValueObservation::try_with_quadratic_penalty(
                SourceId::new("soft-value"),
                soft_location,
                2.0,
                QuadraticPenalty::try_new(4.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    penalty_problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();
    let penalty_fit = penalty_problem.build().unwrap().fit().unwrap();

    let mut statistical_problem = manufactured_problem();
    statistical_problem
        .add(
            FieldValueObservation::try_with_standard_deviation(
                SourceId::new("soft-value"),
                soft_location,
                2.0,
                StandardDeviation::try_new(0.5).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    statistical_problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();
    let statistical_fit = statistical_problem.build().unwrap().fit().unwrap();

    for point in [
        soft_location,
        Point3::try_new(0.0, 0.0, 0.0).unwrap(),
        Point3::try_new(0.6, -0.1, -0.4).unwrap(),
    ] {
        let penalty = penalty_fit.model().evaluate(point).unwrap();
        let statistical = statistical_fit.model().evaluate(point).unwrap();
        assert_close(penalty.value(), statistical.value(), 1.0e-11);
        for (left, right) in penalty
            .gradient()
            .components()
            .into_iter()
            .zip(statistical.gradient().components())
        {
            assert_close(left, right, 1.0e-11);
        }
    }

    let penalty_assessment = &penalty_fit.report().soft_field_values()[0];
    let statistical_assessment = &statistical_fit.report().soft_field_values()[0];
    assert_eq!(
        penalty_assessment.quadratic_penalty(),
        Some(QuadraticPenalty::try_new(4.0).unwrap())
    );
    assert_eq!(penalty_assessment.standard_deviation(), None);
    assert_eq!(statistical_assessment.quadratic_penalty(), None);
    assert_eq!(
        statistical_assessment.standard_deviation(),
        Some(StandardDeviation::try_new(0.5).unwrap())
    );
    assert_close(
        penalty_assessment.loss(),
        statistical_assessment.loss(),
        1.0e-12,
    );
    assert_close(
        penalty_fit.report().total_objective().unwrap(),
        statistical_fit.report().total_objective().unwrap(),
        1.0e-12,
    );
}

#[test]
fn duplicate_soft_evidence_retains_independent_residuals_and_loss_contributions() {
    let location = Point3::try_new(0.2, -0.3, 0.4).unwrap();
    let mut duplicated = manufactured_problem();
    for source in ["soft-a", "soft-b"] {
        duplicated
            .add(
                FieldValueObservation::try_with_quadratic_penalty(
                    SourceId::new(source),
                    location,
                    2.0,
                    QuadraticPenalty::try_new(1.0).unwrap(),
                )
                .unwrap(),
            )
            .unwrap();
    }
    duplicated
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();
    let duplicated_fit = duplicated.build().unwrap().fit().unwrap();

    let mut combined_precision = manufactured_problem();
    combined_precision
        .add(
            FieldValueObservation::try_with_quadratic_penalty(
                SourceId::new("soft-combined"),
                location,
                2.0,
                QuadraticPenalty::try_new(2.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    combined_precision
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();
    let combined_fit = combined_precision.build().unwrap().fit().unwrap();

    let assessments = duplicated_fit.report().soft_field_values();
    assert_eq!(assessments.len(), 2);
    assert_eq!(assessments[0].source_id(), &SourceId::new("soft-a"));
    assert_eq!(assessments[1].source_id(), &SourceId::new("soft-b"));
    assert_close(
        assessments[0].residual(),
        assessments[1].residual(),
        1.0e-13,
    );
    assert_close(assessments[0].loss(), assessments[1].loss(), 1.0e-13);

    let duplicated_sample = duplicated_fit.model().evaluate(location).unwrap();
    let combined_sample = combined_fit.model().evaluate(location).unwrap();
    assert_close(duplicated_sample.value(), combined_sample.value(), 1.0e-11);
    assert_close(
        assessments.iter().map(|assessment| assessment.loss()).sum(),
        combined_fit.report().soft_field_values()[0].loss(),
        1.0e-11,
    );
    assert_close(
        duplicated_fit.report().total_objective().unwrap(),
        combined_fit.report().total_objective().unwrap(),
        1.0e-11,
    );
}

#[test]
fn soft_objective_is_covariant_under_length_and_field_unit_rescaling() {
    const LENGTH_SCALE: f64 = 2.5;
    const FIELD_SCALE: f64 = 4.0;
    const NORMALIZATION: f64 = 3.0;
    const PENALTY: f64 = 2.0;
    let location = [0.2, -0.3, 0.4];

    let mut original = manufactured_problem();
    original
        .add(
            FieldValueObservation::try_with_quadratic_penalty(
                SourceId::new("soft-value"),
                Point3::try_new(location[0], location[1], location[2]).unwrap(),
                2.0,
                QuadraticPenalty::try_new(PENALTY).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    original
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(NORMALIZATION).unwrap())
        .unwrap();
    let original = original.build().unwrap().fit().unwrap();

    let mut rescaled = manufactured_problem_with_scales(LENGTH_SCALE, FIELD_SCALE);
    let [x, y, z] = location.map(|coordinate| coordinate * LENGTH_SCALE);
    rescaled
        .add(
            FieldValueObservation::try_with_quadratic_penalty(
                SourceId::new("soft-value"),
                Point3::try_new(x, y, z).unwrap(),
                FIELD_SCALE * 2.0,
                QuadraticPenalty::try_new(PENALTY / FIELD_SCALE.powi(2)).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    rescaled
        .set_field_energy_normalization(
            FieldEnergyNormalization::try_new(
                NORMALIZATION * LENGTH_SCALE.powi(3) / FIELD_SCALE.powi(2),
            )
            .unwrap(),
        )
        .unwrap();
    let rescaled = rescaled.build().unwrap().fit().unwrap();

    let original_sample = original
        .model()
        .evaluate(Point3::try_new(location[0], location[1], location[2]).unwrap())
        .unwrap();
    let rescaled_sample = rescaled
        .model()
        .evaluate(Point3::try_new(x, y, z).unwrap())
        .unwrap();
    assert_close(
        rescaled_sample.value(),
        FIELD_SCALE * original_sample.value(),
        2.0e-10,
    );
    for (actual, expected) in rescaled_sample
        .gradient()
        .components()
        .into_iter()
        .zip(original_sample.gradient().components())
    {
        assert_close(actual, FIELD_SCALE / LENGTH_SCALE * expected, 2.0e-10);
    }
    assert_close(
        rescaled.report().soft_field_values()[0].residual(),
        FIELD_SCALE * original.report().soft_field_values()[0].residual(),
        2.0e-10,
    );
    assert_close(
        rescaled.report().soft_field_values()[0].loss(),
        original.report().soft_field_values()[0].loss(),
        2.0e-10,
    );
    assert_close(
        rescaled.report().field_energy().unwrap(),
        original.report().field_energy().unwrap(),
        2.0e-10,
    );
    assert_close(
        rescaled.report().total_objective().unwrap(),
        original.report().total_objective().unwrap(),
        2.0e-10,
    );
    for (actual, expected) in rescaled
        .report()
        .hard_relations()
        .iter()
        .zip(original.report().hard_relations())
    {
        assert_eq!(actual.source_id(), expected.source_id());
        assert_close(
            actual.tolerance(),
            FIELD_SCALE * expected.tolerance(),
            2.0e-10,
        );
    }
}
