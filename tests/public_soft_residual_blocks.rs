use georbf::geometry::{
    FieldUnitLabel, GlobalAnisotropyMetric, Handedness, InputCoordinateFrame, LengthUnitLabel,
    Point3, Vector3,
};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{
    CovarianceGroupBuildError, CovarianceGroupBuilder, CovarianceGroupMemberAddError,
    CovarianceMatrix, CovarianceMatrixError, CovarianceResidualDimension, FieldValueObservation,
    GradientObservation, ObservationError, QuadraticPenalty, StandardDeviation,
    TangentDirectionObservation,
};
use georbf::problem::{AddError, BuildError, ProblemInput};
use georbf::relation::AdditiveFieldGauge;
use georbf::{GroupId, ProblemBuilder, SourceId};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).unwrap()
}

fn problem_builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field-unit"),
    )
}

fn assert_requires_normalization<T: ProblemInput>(input: T) {
    let mut builder = problem_builder();
    builder.add(input).unwrap();
    let failure = builder
        .build()
        .expect_err("every soft residual requires an explicit FieldEnergy scale");
    assert_eq!(
        failure.errors(),
        &[BuildError::MissingFieldEnergyNormalization]
    );
    let mut builder = failure.into_builder();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();
    assert!(builder.build().is_ok());
}

fn manufactured_problem() -> ProblemBuilder {
    let mut builder = problem_builder();
    for (index, (support, value)) in [
        ([-1.0, -1.0, -1.0], 0.0),
        ([1.0, -1.0, -1.0], 1.0),
        ([-1.0, 1.0, -1.0], -0.5),
        ([-1.0, -1.0, 1.0], 1.5),
        ([1.0, 1.0, 0.5], 1.625),
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("hard-{index}")),
                    point(support[0], support[1], support[2]),
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn assert_query_components_equivalent(actual: [f64; 3], expected: [f64; 3]) {
    let field_scale = actual
        .into_iter()
        .chain(expected)
        .map(f64::abs)
        .fold(1.0_f64, f64::max);
    for (actual, expected) in actual.into_iter().zip(expected) {
        let sample_reference_scale = actual.abs().max(expected.abs());
        let tolerance = 1.0e-12 * field_scale + 1.0e-11 * sample_reference_scale;
        assert_close(actual, expected, tolerance);
    }
}

#[test]
fn covariance_and_independent_vector_soft_inputs_have_checked_public_boundaries() {
    let covariance =
        CovarianceMatrix::try_new([[0.25, 0.05, 0.0], [0.05, 0.36, -0.02], [0.0, -0.02, 0.49]])
            .expect("the manufactured covariance is finite, symmetric, and SPD");
    assert_eq!(covariance.dimension(), 3);
    assert_eq!(covariance.entry(0, 1), Some(0.05));

    assert_eq!(
        CovarianceMatrix::try_new([[1.0, f64::NAN], [f64::NAN, 1.0]]),
        Err(CovarianceMatrixError::NonFinite { row: 0, column: 1 })
    );
    assert_eq!(
        CovarianceMatrix::try_new([[1.0, 0.25], [0.5, 1.0]]),
        Err(CovarianceMatrixError::NotSymmetric { row: 0, column: 1 })
    );
    assert_eq!(
        CovarianceMatrix::try_new([[1.0, 1.0], [1.0, 1.0]]),
        Err(CovarianceMatrixError::NotPositiveDefinite)
    );

    let location = point(0.25, -0.5, 0.75);
    let gradient = vector(0.5, -1.0, 0.25);
    assert_eq!(
        GradientObservation::try_with_covariance(
            SourceId::new("wrong-dimension"),
            location,
            gradient,
            CovarianceMatrix::try_new([[1.0, 0.0], [0.0, 1.0]]).unwrap(),
        ),
        Err(ObservationError::CovarianceDimensionMismatch {
            expected: 3,
            actual: 2,
        })
    );

    assert_requires_normalization(GradientObservation::with_quadratic_penalty(
        SourceId::new("soft-gradient-penalty"),
        location,
        gradient,
        QuadraticPenalty::try_new(2.0).unwrap(),
    ));
    assert_requires_normalization(GradientObservation::with_standard_deviation(
        SourceId::new("soft-gradient-isotropic"),
        location,
        gradient,
        StandardDeviation::try_new(0.3).unwrap(),
    ));
    assert_requires_normalization(
        GradientObservation::try_with_covariance(
            SourceId::new("soft-gradient-covariance"),
            location,
            gradient,
            covariance,
        )
        .unwrap(),
    );
    assert_requires_normalization(
        TangentDirectionObservation::try_with_quadratic_penalty(
            SourceId::new("soft-tangent-penalty"),
            location,
            vector(1.0, 1.0, 0.0),
            QuadraticPenalty::try_new(3.0).unwrap(),
        )
        .unwrap(),
    );
    assert_requires_normalization(
        TangentDirectionObservation::try_with_standard_deviation(
            SourceId::new("soft-tangent-statistical"),
            location,
            vector(1.0, -1.0, 0.0),
            StandardDeviation::try_new(0.2).unwrap(),
        )
        .unwrap(),
    );
}

#[test]
fn euclidean_soft_gradient_is_one_ordered_vector_residual_block() {
    let location = point(0.2, -0.3, 0.4);
    let target = vector(1.5, -0.75, 0.25);
    let penalty = QuadraticPenalty::try_new(2.0).unwrap();
    let mut builder = manufactured_problem();
    builder
        .set_global_anisotropy_metric(
            GlobalAnisotropyMetric::try_from_matrix([
                [2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.5],
            ])
            .unwrap(),
        )
        .unwrap();
    builder
        .add(GradientObservation::with_quadratic_penalty(
            SourceId::new("soft-gradient"),
            location,
            target,
            penalty,
        ))
        .unwrap();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();

    let fit = builder
        .build()
        .unwrap()
        .fit()
        .expect("a vector quadratic residual stays on the symmetric KKT path");
    let report = fit.report();
    assert_eq!(report.problem_size().scalar_soft_relations(), 3);
    assert_eq!(report.problem_size().quadratic_objective_terms(), 1);
    assert!(report.direct_input_conflicts().is_empty());

    let assessment = &report.soft_gradients()[0];
    assert_eq!(assessment.source_id(), &SourceId::new("soft-gradient"));
    assert_eq!(assessment.target(), target);
    assert_eq!(assessment.quadratic_penalty(), Some(penalty));
    assert_eq!(assessment.covariance(), None);
    let sampled = fit.model().evaluate(location).unwrap().gradient();
    assert_query_components_equivalent(
        assessment.recovered_gradient().components(),
        sampled.components(),
    );
    let residual = assessment.residual().components();
    for ((actual, recovered), expected) in residual
        .into_iter()
        .zip(sampled.components())
        .zip(target.components())
    {
        assert_close(actual, recovered - expected, 1.0e-13);
    }
    let expected_loss = 0.5
        * penalty.weight()
        * residual
            .into_iter()
            .map(|component| component * component)
            .sum::<f64>();
    assert_close(assessment.loss(), expected_loss, 1.0e-12);
    assert_close(
        report.total_objective().unwrap(),
        0.5 * report.field_energy().unwrap() + assessment.loss(),
        1.0e-11,
    );
}

#[test]
fn gradient_covariance_uses_cross_terms_and_recovers_the_physical_vector() {
    let location = point(0.2, -0.3, 0.4);
    let target = vector(1.5, -0.75, 0.25);
    let covariance =
        CovarianceMatrix::try_new([[1.0, 0.5, 0.0], [0.5, 1.0, 0.0], [0.0, 0.0, 4.0]]).unwrap();
    let mut builder = manufactured_problem();
    builder
        .add(
            GradientObservation::try_with_covariance(
                SourceId::new("statistical-gradient"),
                location,
                target,
                covariance.clone(),
            )
            .unwrap(),
        )
        .unwrap();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();

    let fit = builder.build().unwrap().fit().unwrap();
    let assessment = &fit.report().soft_gradients()[0];
    assert_eq!(assessment.quadratic_penalty(), None);
    assert_eq!(assessment.covariance(), Some(&covariance));
    assert!(assessment.whitening_round_trip_error() <= 1.0e-11);
    assert!(
        fit.report()
            .canonical_acceptance()
            .unwrap()
            .whitening_round_trip_error()
            .unwrap()
            <= 1.0e-11
    );
    let [rx, ry, rz] = assessment.residual().components();
    let expected_loss = 0.5
        * (rx * (4.0 / 3.0 * rx - 2.0 / 3.0 * ry)
            + ry * (-2.0 / 3.0 * rx + 4.0 / 3.0 * ry)
            + rz * 0.25 * rz);
    assert_close(assessment.loss(), expected_loss, 1.0e-11);
    assert_close(
        assessment
            .whitened_residual()
            .iter()
            .map(|component| component * component)
            .sum::<f64>(),
        2.0 * assessment.loss(),
        1.0e-11,
    );
}

#[test]
fn soft_tangent_remains_one_zero_directional_derivative_residual() {
    let location = point(0.2, -0.3, 0.4);
    let direction = vector(2.0, -1.0, 2.0);
    let standard_deviation = StandardDeviation::try_new(0.4).unwrap();
    let mut builder = manufactured_problem();
    builder
        .add(
            TangentDirectionObservation::try_with_standard_deviation(
                SourceId::new("soft-tangent"),
                location,
                direction,
                standard_deviation,
            )
            .unwrap(),
        )
        .unwrap();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();

    let fit = builder.build().unwrap().fit().unwrap();
    let report = fit.report();
    assert_eq!(report.problem_size().scalar_soft_relations(), 1);
    assert_eq!(report.problem_size().quadratic_objective_terms(), 1);
    let assessment = &report.soft_tangents()[0];
    assert_eq!(assessment.source_id(), &SourceId::new("soft-tangent"));
    assert_eq!(assessment.target(), 0.0);
    assert_eq!(assessment.quadratic_penalty(), None);
    assert_eq!(assessment.standard_deviation(), Some(standard_deviation));

    let gradient = fit
        .model()
        .evaluate(location)
        .unwrap()
        .gradient()
        .components();
    let unit = [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0];
    let directional_derivative = gradient
        .into_iter()
        .zip(unit)
        .map(|(component, direction)| component * direction)
        .sum::<f64>();
    assert_close(
        assessment.recovered_directional_derivative(),
        directional_derivative,
        1.0e-12,
    );
    assert_close(assessment.residual(), directional_derivative, 1.0e-12);
    assert_close(
        assessment.loss(),
        0.5 * (directional_derivative / standard_deviation.value()).powi(2),
        1.0e-12,
    );
}

#[test]
fn zero_gradient_satisfies_a_soft_tangent_without_polarity_or_slope() {
    let mut builder = problem_builder();
    for (index, support) in [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, 1.0, 0.5],
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("constant-{index}")),
                    point(support[0], support[1], support[2]),
                    2.0,
                )
                .unwrap(),
            )
            .unwrap();
    }
    let location = point(0.2, -0.3, 0.4);
    builder
        .add(
            TangentDirectionObservation::try_with_quadratic_penalty(
                SourceId::new("zero-gradient-tangent"),
                location,
                vector(1.0, 2.0, 3.0),
                QuadraticPenalty::try_new(2.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let fit = builder.build().unwrap().fit().unwrap();
    let sample = fit.model().evaluate(location).unwrap();
    for component in sample.gradient().components() {
        assert_close(component, 0.0, 1.0e-10);
    }
    let assessment = &fit.report().soft_tangents()[0];
    assert_close(assessment.residual(), 0.0, 1.0e-10);
    assert_close(assessment.loss(), 0.0, 1.0e-10);
}

#[test]
fn covariance_groups_are_complete_same_dimension_atomic_problem_inputs() {
    let empty = CovarianceGroupBuilder::new(GroupId::new("empty"));
    let empty_failure = empty
        .build(CovarianceMatrix::try_new([[1.0]]).unwrap())
        .expect_err("an empty group is incomplete");
    assert_eq!(
        empty_failure.error(),
        &CovarianceGroupBuildError::EmptyGroup
    );
    let (mut repaired_empty, original_covariance) = empty_failure.into_parts();
    repaired_empty
        .add_field_value_member(SourceId::new("repaired/member"), point(0.0, 0.0, 0.0), 1.0)
        .unwrap();
    assert!(repaired_empty.build(original_covariance).is_ok());

    let mut draft = CovarianceGroupBuilder::new(GroupId::new("derivative-group"));
    draft
        .add_gradient_member(
            SourceId::new("gradient-z"),
            point(0.2, -0.3, 0.4),
            vector(1.5, -0.75, 0.25),
        )
        .unwrap();
    assert_eq!(
        draft.add_field_value_member(SourceId::new("wrong-dimension"), point(0.0, 0.0, 0.0), 1.0,),
        Err(CovarianceGroupMemberAddError::DimensionMismatch {
            expected: CovarianceResidualDimension::FieldValuePerLength,
            actual: CovarianceResidualDimension::FieldValue,
        })
    );
    assert_eq!(
        draft.add_tangent_member(
            SourceId::new("gradient-z"),
            point(0.0, 0.0, 0.0),
            vector(1.0, 0.0, 0.0),
        ),
        Err(CovarianceGroupMemberAddError::DuplicateSourceId {
            source_id: SourceId::new("gradient-z"),
        })
    );
    draft
        .add_tangent_member(
            SourceId::new("tangent-a"),
            point(-0.1, 0.2, -0.4),
            vector(1.0, 1.0, 0.0),
        )
        .expect("both rejected mutations left the draft repairable");

    let mut wrong_dimension = CovarianceGroupBuilder::new(GroupId::new("wrong-covariance"));
    wrong_dimension
        .add_gradient_member(
            SourceId::new("wrong-covariance/member"),
            point(0.0, 0.0, 0.0),
            vector(0.0, 0.0, 0.0),
        )
        .unwrap();
    let wrong_dimension_failure = wrong_dimension
        .build(CovarianceMatrix::try_new([[1.0, 0.0], [0.0, 1.0]]).unwrap())
        .expect_err("the covariance does not cover all three gradient components");
    assert_eq!(
        wrong_dimension_failure.error(),
        &CovarianceGroupBuildError::CovarianceDimensionMismatch {
            expected: 3,
            actual: 2,
        }
    );
    let (wrong_dimension, rejected_covariance) = wrong_dimension_failure.into_parts();
    assert_eq!(rejected_covariance.dimension(), 2);
    assert!(
        wrong_dimension
            .build(
                CovarianceMatrix::try_new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],])
                    .unwrap()
            )
            .is_ok(),
        "a failed build retains the complete draft for repair"
    );

    let covariance = CovarianceMatrix::try_new([
        [1.0, 0.0, 0.0, 0.25],
        [0.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 3.0, 0.0],
        [0.25, 0.0, 0.0, 1.0],
    ])
    .unwrap();
    let group = draft.build(covariance).unwrap();
    assert_eq!(group.group_id(), &GroupId::new("derivative-group"));
    assert_eq!(group.member_count(), 2);
    assert_eq!(group.scalar_residual_count(), 4);

    let mut builder = manufactured_problem();
    builder.add(group).unwrap();

    let mut duplicate_source = CovarianceGroupBuilder::new(GroupId::new("new-group"));
    duplicate_source
        .add_tangent_member(
            SourceId::new("gradient-z"),
            point(0.0, 0.0, 0.0),
            vector(1.0, 0.0, 0.0),
        )
        .unwrap();
    assert_eq!(
        builder.add(
            duplicate_source
                .build(CovarianceMatrix::try_new([[1.0]]).unwrap())
                .unwrap()
        ),
        Err(AddError::DuplicateSourceId {
            source_id: SourceId::new("gradient-z"),
        })
    );

    let mut duplicate_group = CovarianceGroupBuilder::new(GroupId::new("derivative-group"));
    duplicate_group
        .add_tangent_member(
            SourceId::new("new-source"),
            point(0.0, 0.0, 0.0),
            vector(1.0, 0.0, 0.0),
        )
        .unwrap();
    assert_eq!(
        builder.add(
            duplicate_group
                .build(CovarianceMatrix::try_new([[1.0]]).unwrap())
                .unwrap()
        ),
        Err(AddError::DuplicateGroupId {
            group_id: GroupId::new("derivative-group"),
        })
    );

    let failure = builder
        .build()
        .expect_err("the accepted covariance group is a soft relation");
    assert_eq!(
        failure.errors(),
        &[BuildError::MissingFieldEnergyNormalization]
    );
    let mut builder = failure.into_builder();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();
    let snapshot = builder.build().unwrap();
    assert_eq!(snapshot.covariance_group_count(), 1);
    assert_eq!(snapshot.source_count(), 7);
}

#[test]
fn cross_member_covariance_reports_one_group_objective_without_member_losses() {
    let gradient_location = point(0.2, -0.3, 0.4);
    let tangent_location = point(-0.1, 0.2, -0.4);
    let mut draft = CovarianceGroupBuilder::new(GroupId::new("derivative-group"));
    draft
        .add_gradient_member(
            SourceId::new("gradient-z"),
            gradient_location,
            vector(1.5, -0.75, 0.25),
        )
        .unwrap();
    draft
        .add_tangent_member(
            SourceId::new("tangent-a"),
            tangent_location,
            vector(1.0, 1.0, 0.0),
        )
        .unwrap();
    let covariance = CovarianceMatrix::try_new([
        [1.0, 0.0, 0.0, 0.5],
        [0.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 3.0, 0.0],
        [0.5, 0.0, 0.0, 1.0],
    ])
    .unwrap();

    let mut builder = manufactured_problem();
    builder
        .add(draft.build(covariance.clone()).unwrap())
        .unwrap();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();
    let fit = builder.build().unwrap().fit().unwrap();
    let report = fit.report();
    assert_eq!(report.problem_size().scalar_soft_relations(), 4);
    assert_eq!(report.problem_size().quadratic_objective_terms(), 1);
    assert!(report.soft_gradients().is_empty());
    assert!(report.soft_tangents().is_empty());

    let group = &report.covariance_groups()[0];
    assert_eq!(group.group_id(), &GroupId::new("derivative-group"));
    assert_eq!(group.covariance(), &covariance);
    assert_eq!(
        group
            .members()
            .iter()
            .map(|member| member.source_id().as_str())
            .collect::<Vec<_>>(),
        ["gradient-z", "tangent-a"]
    );
    assert_eq!(group.members()[0].residual_components().len(), 3);
    assert_eq!(group.members()[1].residual_components().len(), 1);
    assert!(group.whitening_round_trip_error() <= 1.0e-11);

    let residual = group
        .members()
        .iter()
        .flat_map(|member| member.residual_components().iter().copied())
        .collect::<Vec<_>>();
    let [rx, ry, rz, rt] = <[f64; 4]>::try_from(residual).unwrap();
    let expected_loss = 0.5
        * (rx * (4.0 / 3.0 * rx - 2.0 / 3.0 * rt)
            + 0.5 * ry * ry
            + (1.0 / 3.0) * rz * rz
            + rt * (-2.0 / 3.0 * rx + 4.0 / 3.0 * rt));
    assert_close(group.objective_contribution(), expected_loss, 1.0e-11);
    assert_close(
        report.total_objective().unwrap(),
        0.5 * report.field_energy().unwrap() + group.objective_contribution(),
        1.0e-11,
    );

    let gradient = fit
        .model()
        .evaluate(gradient_location)
        .unwrap()
        .gradient()
        .components();
    assert_query_components_equivalent(
        <[f64; 3]>::try_from(group.members()[0].recovered_components()).unwrap(),
        gradient,
    );
    let tangent_gradient = fit
        .model()
        .evaluate(tangent_location)
        .unwrap()
        .gradient()
        .components();
    let recovered_tangent = (tangent_gradient[0] + tangent_gradient[1]) / 2.0_f64.sqrt();
    assert_close(
        group.members()[1].recovered_components()[0],
        recovered_tangent,
        1.0e-12,
    );
}

#[test]
fn conflicting_soft_gradients_increase_objective_without_a_conflict_diagnosis() {
    let location = point(0.2, -0.3, 0.4);
    let mut builder = manufactured_problem();
    for (source, target) in [
        ("soft-conflict-a", vector(2.0, 0.0, 0.0)),
        ("soft-conflict-b", vector(-2.0, 0.0, 0.0)),
    ] {
        builder
            .add(GradientObservation::with_quadratic_penalty(
                SourceId::new(source),
                location,
                target,
                QuadraticPenalty::try_new(1.0).unwrap(),
            ))
            .unwrap();
    }
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(3.0).unwrap())
        .unwrap();

    let fit = builder
        .build()
        .unwrap()
        .fit()
        .expect("soft contradictions are objective tradeoffs, not infeasibility");
    assert!(fit.report().direct_input_conflicts().is_empty());
    assert_eq!(fit.report().soft_gradients().len(), 2);
    assert!(
        fit.report()
            .soft_gradients()
            .iter()
            .map(|assessment| assessment.loss())
            .sum::<f64>()
            > 0.0
    );
}

#[test]
fn covariance_group_field_value_anchor_keeps_an_additive_gauge_out_of_the_solver() {
    let anchor_location = point(0.0, 0.0, 0.0);
    let mut group = CovarianceGroupBuilder::new(GroupId::new("absolute-statistical-group"));
    group
        .add_field_value_member(SourceId::new("soft-anchor"), anchor_location, 2.0)
        .unwrap();

    let mut builder = problem_builder();
    builder
        .add(
            group
                .build(CovarianceMatrix::try_new([[1.0]]).unwrap())
                .unwrap(),
        )
        .unwrap();
    builder
        .add(GradientObservation::new(
            SourceId::new("complete-gradient"),
            anchor_location,
            vector(1.0, -0.5, 0.25),
        ))
        .unwrap();
    builder
        .add(
            AdditiveFieldGauge::at_point(SourceId::new("reporting-gauge"), anchor_location, 2.0)
                .unwrap(),
        )
        .unwrap();
    builder
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let fit = builder.build().unwrap().fit().unwrap();
    let gauge = fit
        .report()
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id() == &SourceId::new("reporting-gauge"))
        .expect("the report retains the verification-only gauge");
    assert!(
        gauge.scaled_kkt_tolerance().is_none(),
        "a statistical field anchor already fixes the additive representative"
    );
}

#[derive(Clone, Copy)]
enum CovariantFrameCase {
    Original,
    Rotated,
    Reflected,
}

fn orthogonal_matrix(case: CovariantFrameCase) -> [[f64; 3]; 3] {
    match case {
        CovariantFrameCase::Original => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        CovariantFrameCase::Rotated => [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        CovariantFrameCase::Reflected => [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    }
}

fn matrix_vector(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|row| {
        (0..3)
            .map(|column| matrix[row][column] * vector[column])
            .sum()
    })
}

fn transform_covariance(
    matrix: [[f64; 3]; 3],
    covariance: [[f64; 3]; 3],
    factor: f64,
) -> [[f64; 3]; 3] {
    std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            factor
                * (0..3)
                    .flat_map(|left| {
                        (0..3).map(move |right| {
                            matrix[row][left] * covariance[left][right] * matrix[column][right]
                        })
                    })
                    .sum::<f64>()
        })
    })
}

fn covariant_soft_gradient_problem(case: CovariantFrameCase) -> georbf::ProblemSnapshot {
    const LENGTH_SCALE: f64 = 2.5;
    const FIELD_SCALE: f64 = 4.0;
    let transformed = !matches!(case, CovariantFrameCase::Original);
    let matrix = orthogonal_matrix(case);
    let (scale, field_scale) = if transformed {
        (LENGTH_SCALE, FIELD_SCALE)
    } else {
        (1.0, 1.0)
    };
    let frame = InputCoordinateFrame::try_new(
        match case {
            CovariantFrameCase::Original => ["east", "north", "elevation"],
            CovariantFrameCase::Rotated => ["rotated-north", "rotated-east", "rotated-up"],
            CovariantFrameCase::Reflected => ["reflected-north", "reflected-east", "reflected-up"],
        },
        if matches!(case, CovariantFrameCase::Reflected) {
            Handedness::Left
        } else {
            Handedness::Right
        },
        LengthUnitLabel::new(if transformed { "scaled-m" } else { "m" }),
    )
    .unwrap();
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("field-unit"));
    let transform_point = |components: [f64; 3]| {
        if transformed {
            let rotated = matrix_vector(matrix, components);
            [
                scale * rotated[0] + 10.0,
                scale * rotated[1] - 3.0,
                scale * rotated[2] + 4.0,
            ]
        } else {
            components
        }
    };
    let transform_gradient_components = |components: [f64; 3]| {
        matrix_vector(matrix, components).map(|component| field_scale / scale * component)
    };
    for (index, (support, value)) in [
        ([-1.0, -1.0, -1.0], 0.0),
        ([1.0, -1.0, -1.0], 1.0),
        ([-1.0, 1.0, -1.0], -0.5),
        ([-1.0, -1.0, 1.0], 1.5),
        ([1.0, 1.0, 0.5], 1.625),
    ]
    .into_iter()
    .enumerate()
    {
        let support = transform_point(support);
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("hard-{index}")),
                    point(support[0], support[1], support[2]),
                    field_scale * value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    let original_covariance = [[1.0, 0.2, 0.1], [0.2, 2.0, -0.15], [0.1, -0.15, 3.0]];
    let covariance =
        transform_covariance(matrix, original_covariance, (field_scale / scale).powi(2));
    let location = transform_point([0.2, -0.3, 0.4]);
    let target = transform_gradient_components([1.5, -0.75, 0.25]);
    builder
        .add(
            GradientObservation::try_with_covariance(
                SourceId::new("statistical-gradient"),
                point(location[0], location[1], location[2]),
                vector(target[0], target[1], target[2]),
                CovarianceMatrix::try_new(covariance).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    let original_metric = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]];
    builder
        .set_global_anisotropy_metric(
            GlobalAnisotropyMetric::try_from_matrix(transform_covariance(
                matrix,
                original_metric,
                1.0,
            ))
            .unwrap(),
        )
        .unwrap();
    builder
        .set_field_energy_normalization(
            FieldEnergyNormalization::try_new(3.0 * scale.powi(3) / field_scale.powi(2)).unwrap(),
        )
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn soft_covariance_objective_is_covariant_under_rotation_reflection_and_unit_scaling() {
    let original = covariant_soft_gradient_problem(CovariantFrameCase::Original)
        .fit()
        .unwrap();
    for case in [CovariantFrameCase::Rotated, CovariantFrameCase::Reflected] {
        let transformed = covariant_soft_gradient_problem(case).fit().unwrap();
        let matrix = orthogonal_matrix(case);
        let query = [0.35, -0.2, 0.1];
        let rotated = matrix_vector(matrix, query);
        let transformed_query = [
            2.5 * rotated[0] + 10.0,
            2.5 * rotated[1] - 3.0,
            2.5 * rotated[2] + 4.0,
        ];
        let original_sample = original
            .model()
            .evaluate(point(query[0], query[1], query[2]))
            .unwrap();
        let transformed_sample = transformed
            .model()
            .evaluate(point(
                transformed_query[0],
                transformed_query[1],
                transformed_query[2],
            ))
            .unwrap();
        assert_close(
            transformed_sample.value(),
            4.0 * original_sample.value(),
            2.0e-8,
        );
        let expected_gradient = matrix_vector(matrix, original_sample.gradient().components())
            .map(|component| 4.0 / 2.5 * component);
        for (actual, expected) in transformed_sample
            .gradient()
            .components()
            .into_iter()
            .zip(expected_gradient)
        {
            assert_close(actual, expected, 2.0e-8);
        }

        let original_assessment = &original.report().soft_gradients()[0];
        let transformed_assessment = &transformed.report().soft_gradients()[0];
        let expected_residual = matrix_vector(matrix, original_assessment.residual().components())
            .map(|component| 4.0 / 2.5 * component);
        for (actual, expected) in transformed_assessment
            .residual()
            .components()
            .into_iter()
            .zip(expected_residual)
        {
            assert_close(actual, expected, 2.0e-8);
        }
        assert_close(
            transformed_assessment.loss(),
            original_assessment.loss(),
            2.0e-8,
        );
        assert_close(
            transformed.report().field_energy().unwrap(),
            original.report().field_energy().unwrap(),
            2.0e-8,
        );
        assert_close(
            transformed.report().total_objective().unwrap(),
            original.report().total_objective().unwrap(),
            2.0e-8,
        );
    }
}
