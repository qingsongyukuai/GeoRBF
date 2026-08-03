use georbf::diagnostics::ProblemDiagnosis;
use georbf::fit::BoundActiveState;
use georbf::geometry::{
    FieldUnitLabel, GlobalAnisotropyMetric, Handedness, InputCoordinateFrame, LengthUnitLabel,
    Point3, Vector3,
};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{
    AxialNormalObservation, DirectedNormalObservation, GradientObservation, MinimumNormalSlope,
    MinimumNormalSlopeEnforcement, MinimumNormalSlopeError, NormalDirectionEnforcement,
    NormalObservationError, QuadraticPenalty, StandardDeviation,
};
use georbf::problem::BuildError;
use georbf::relation::SharedLevelSetBuilder;
use georbf::relation::{LinearViolationPenalty, PolarityResolution, PolaritySelection};
use georbf::{GroupId, ProblemBuilder, SourceId};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).unwrap()
}

fn builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("field"),
    )
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 5.0e-7 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn add_affine_values(problem: &mut ProblemBuilder, gradient: [f64; 3]) {
    for (index, coordinates) in [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, 1.0, 0.5],
    ]
    .into_iter()
    .enumerate()
    {
        let [x, y, z] = coordinates;
        let value = 4.0 + gradient[0] * x + gradient[1] * y + gradient[2] * z;
        problem
            .add(
                georbf::observation::FieldValueObservation::try_new(
                    SourceId::new(format!("value-{index}")),
                    point(x, y, z),
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }
}

fn add_affine_shared_level(problem: &mut ProblemBuilder) {
    let mut level = SharedLevelSetBuilder::new(GroupId::new("normal-level"));
    level
        .add_member(SourceId::new("normal-level-a"), point(0.0, 0.0, 0.0))
        .unwrap();
    level
        .add_member(SourceId::new("normal-level-b"), point(2.0, -1.0, 0.0))
        .unwrap();
    problem.add(level.build().unwrap()).unwrap();
}

#[test]
fn checked_normal_inputs_preserve_directed_polarity_and_axial_identity() {
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let expected = if value.is_finite() {
            MinimumNormalSlopeError::NotPositive
        } else {
            MinimumNormalSlopeError::NotFinite
        };
        assert_eq!(MinimumNormalSlope::try_new(value), Err(expected));
    }
    let slope = MinimumNormalSlope::try_new(0.25).unwrap();
    assert_eq!(slope.value(), 0.25);

    let directed = DirectedNormalObservation::try_new(
        SourceId::new("directed"),
        point(1.0, 2.0, 3.0),
        vector(2.0, -3.0, 6.0),
        slope,
    )
    .unwrap();
    assert_eq!(
        directed.direction().components(),
        [2.0 / 7.0, -3.0 / 7.0, 6.0 / 7.0]
    );
    assert!(!directed.direction_enforcement().is_soft());
    assert!(!directed.minimum_slope_enforcement().is_soft());

    let scaled = DirectedNormalObservation::try_new(
        SourceId::new("scaled"),
        point(1.0, 2.0, 3.0),
        vector(20.0, -30.0, 60.0),
        slope,
    )
    .unwrap();
    assert_eq!(scaled.direction(), directed.direction());
    let opposite = DirectedNormalObservation::try_new(
        SourceId::new("opposite"),
        point(1.0, 2.0, 3.0),
        vector(-2.0, 3.0, -6.0),
        slope,
    )
    .unwrap();
    assert_eq!(
        opposite.direction().components(),
        [-2.0 / 7.0, 3.0 / 7.0, -6.0 / 7.0]
    );

    let axial = AxialNormalObservation::try_new(
        SourceId::new("axial"),
        point(1.0, 2.0, 3.0),
        vector(2.0, -3.0, 6.0),
        slope,
    )
    .unwrap();
    let reversed_axial = AxialNormalObservation::try_new(
        SourceId::new("axial-reversed"),
        point(1.0, 2.0, 3.0),
        vector(-20.0, 30.0, -60.0),
        slope,
    )
    .unwrap();
    assert_eq!(axial.axis(), reversed_axial.axis());
    assert_eq!(axial.input_axis(), directed.direction());
    assert_eq!(reversed_axial.input_axis(), opposite.direction());

    let zero_direction = vector(0.0, -0.0, 0.0);
    assert_eq!(
        DirectedNormalObservation::try_new(
            SourceId::new("zero-directed"),
            point(0.0, 0.0, 0.0),
            zero_direction,
            slope,
        ),
        Err(NormalObservationError::ZeroDirection)
    );
    assert_eq!(
        AxialNormalObservation::try_new(
            SourceId::new("zero-axial"),
            point(0.0, 0.0, 0.0),
            zero_direction,
            slope,
        ),
        Err(NormalObservationError::ZeroDirection)
    );

    for direction in [
        vector(f64::MAX, -f64::MAX, 0.0),
        vector(f64::MIN_POSITIVE, f64::MIN_POSITIVE, 0.0),
        vector(f64::from_bits(1), 0.0, 0.0),
    ] {
        let unit = DirectedNormalObservation::try_new(
            SourceId::new("extreme"),
            point(0.0, 0.0, 0.0),
            direction,
            slope,
        )
        .unwrap()
        .direction()
        .components();
        assert!(unit.iter().all(|component| component.is_finite()));
        assert!(
            (unit
                .into_iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt()
                - 1.0)
                .abs()
                <= 1.0e-15
        );
    }
}

#[test]
fn normal_direction_and_slope_enforcement_are_independently_typed() {
    let direction_penalty = QuadraticPenalty::try_new(2.0).unwrap();
    let slope_penalty = LinearViolationPenalty::try_new(3.0).unwrap();
    let normal = DirectedNormalObservation::try_with_enforcement(
        SourceId::new("soft-normal"),
        point(0.0, 0.0, 0.0),
        vector(0.0, 0.0, 4.0),
        NormalDirectionEnforcement::with_quadratic_penalty(direction_penalty),
        MinimumNormalSlope::try_new(1.5).unwrap(),
        MinimumNormalSlopeEnforcement::with_linear_violation_penalty(slope_penalty),
    )
    .unwrap();
    assert_eq!(normal.direction().components(), [0.0, 0.0, 1.0]);
    assert_eq!(
        normal.direction_enforcement().quadratic_penalty(),
        Some(direction_penalty)
    );
    assert_eq!(normal.direction_enforcement().standard_deviation(), None);
    assert_eq!(
        normal
            .minimum_slope_enforcement()
            .linear_violation_penalty(),
        Some(slope_penalty)
    );
    assert_eq!(normal.minimum_slope_enforcement().quadratic_penalty(), None);

    let statistical = NormalDirectionEnforcement::with_standard_deviation(
        StandardDeviation::try_new(0.2).unwrap(),
    );
    assert!(statistical.is_soft());
    assert_eq!(statistical.standard_deviation().unwrap().value(), 0.2);
}

#[test]
fn polarity_resolution_is_an_independent_forward_referenceable_input() {
    let resolution = PolarityResolution::new(
        SourceId::new("field-decision-17"),
        SourceId::new("axial-survey-9"),
        PolaritySelection::AgainstInputAxis,
    );
    assert_eq!(resolution.source_id(), &SourceId::new("field-decision-17"));
    assert_eq!(
        resolution.axial_normal_source_id(),
        &SourceId::new("axial-survey-9")
    );
    assert_eq!(resolution.selection(), PolaritySelection::AgainstInputAxis);

    let mut problem = builder();
    problem.add(resolution).unwrap();
}

#[test]
fn unresolved_axial_normal_is_preserved_and_rejected_before_backend_execution() {
    let slope = MinimumNormalSlope::try_new(0.5).unwrap();
    let axial = AxialNormalObservation::try_new(
        SourceId::new("axial-unresolved"),
        point(1.0, -2.0, 3.0),
        vector(-2.0, 3.0, -6.0),
        slope,
    )
    .unwrap();
    let mut problem = builder();
    problem.add(axial.clone()).unwrap();
    let snapshot = problem.build().unwrap();
    assert_eq!(snapshot.axial_normal_count(), 1);
    assert_eq!(snapshot.polarity_resolution_count(), 0);
    assert_eq!(
        snapshot
            .axial_normal_observation(&SourceId::new("axial-unresolved"))
            .unwrap()
            .input_axis(),
        axial.input_axis()
    );

    let failure = snapshot.fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::UnresolvedSemantics);
    let evidence = failure.report().unresolved_axial_normals();
    assert_eq!(evidence.len(), 1);
    assert_eq!(evidence[0].source_id(), axial.source_id());
    assert_eq!(evidence[0].input_axis(), axial.input_axis());
    assert!(!evidence[0].backend_invoked());
    assert!(failure.report().attempts().is_empty());
}

#[test]
fn polarity_resolution_allows_forward_reference_and_rejects_dangling_or_repeated_decisions() {
    let axial_source = SourceId::new("axial-forward");
    let first_source = SourceId::new("resolution-first");
    let mut problem = builder();
    problem
        .add(PolarityResolution::new(
            first_source.clone(),
            axial_source.clone(),
            PolaritySelection::AlongInputAxis,
        ))
        .unwrap();
    problem
        .add(
            AxialNormalObservation::try_new(
                axial_source.clone(),
                point(0.0, 0.0, 0.0),
                vector(0.0, 0.0, 2.0),
                MinimumNormalSlope::try_new(0.25).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    let snapshot = problem.build().unwrap();
    assert_eq!(snapshot.axial_normal_count(), 1);
    assert_eq!(snapshot.polarity_resolution_count(), 1);
    assert_eq!(
        snapshot
            .polarity_resolution(&first_source)
            .unwrap()
            .selection(),
        PolaritySelection::AlongInputAxis
    );

    let duplicate_source = SourceId::new("resolution-duplicate");
    let conflict_source = SourceId::new("resolution-conflict");
    let mut repeated = builder();
    repeated
        .add(
            AxialNormalObservation::try_new(
                axial_source.clone(),
                point(0.0, 0.0, 0.0),
                vector(0.0, 0.0, 1.0),
                MinimumNormalSlope::try_new(0.25).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    for resolution in [
        PolarityResolution::new(
            first_source.clone(),
            axial_source.clone(),
            PolaritySelection::AlongInputAxis,
        ),
        PolarityResolution::new(
            duplicate_source.clone(),
            axial_source.clone(),
            PolaritySelection::AlongInputAxis,
        ),
        PolarityResolution::new(
            conflict_source.clone(),
            axial_source.clone(),
            PolaritySelection::AgainstInputAxis,
        ),
    ] {
        repeated.add(resolution).unwrap();
    }
    let repeated_failure = repeated.build().unwrap_err();
    assert_eq!(
        repeated_failure.errors(),
        &[
            BuildError::DuplicatePolarityResolution {
                axial_normal_source_id: axial_source.clone(),
                selection: PolaritySelection::AlongInputAxis,
                resolution_source_ids: vec![duplicate_source.clone(), first_source.clone()],
            },
            BuildError::ConflictingPolarityResolution {
                axial_normal_source_id: axial_source.clone(),
                along_resolution_source_ids: vec![duplicate_source, first_source],
                against_resolution_source_ids: vec![conflict_source],
            },
        ]
    );
    let expected_repeated_errors = repeated_failure.errors().to_vec();
    let mut reversed = builder();
    reversed
        .add(
            AxialNormalObservation::try_new(
                axial_source.clone(),
                point(0.0, 0.0, 0.0),
                vector(0.0, 0.0, 1.0),
                MinimumNormalSlope::try_new(0.25).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    for resolution in [
        PolarityResolution::new(
            SourceId::new("resolution-conflict"),
            axial_source.clone(),
            PolaritySelection::AgainstInputAxis,
        ),
        PolarityResolution::new(
            SourceId::new("resolution-duplicate"),
            axial_source.clone(),
            PolaritySelection::AlongInputAxis,
        ),
        PolarityResolution::new(
            SourceId::new("resolution-first"),
            axial_source.clone(),
            PolaritySelection::AlongInputAxis,
        ),
    ] {
        reversed.add(resolution).unwrap();
    }
    assert_eq!(
        reversed.build().unwrap_err().errors(),
        expected_repeated_errors
    );

    let missing_axial_source = SourceId::new("missing-axial");
    let dangling_source = SourceId::new("dangling-resolution");
    let mut dangling = builder();
    dangling
        .add(PolarityResolution::new(
            dangling_source.clone(),
            missing_axial_source.clone(),
            PolaritySelection::AlongInputAxis,
        ))
        .unwrap();
    let failure = dangling.build().unwrap_err();
    assert_eq!(
        failure.errors(),
        &[BuildError::UnknownAxialNormalReference {
            resolution_source_id: dangling_source,
            axial_normal_source_id: missing_axial_source,
        }]
    );
}

#[test]
fn hard_directed_normal_reports_projection_and_independent_slope_relation() {
    let mut problem = builder();
    add_affine_values(&mut problem, [1.0, 2.0, 2.0]);
    problem
        .add(
            DirectedNormalObservation::try_new(
                SourceId::new("hard-normal"),
                point(0.25, -0.5, 0.75),
                vector(1.0, 2.0, 2.0),
                MinimumNormalSlope::try_new(2.5).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();

    let fit = problem.build().unwrap().fit().unwrap();
    let assessment = &fit.report().directed_normals()[0];
    assert_eq!(assessment.source_id(), &SourceId::new("hard-normal"));
    assert_eq!(
        assessment.direction_semantic_role().as_str(),
        "directed-normal/direction-projection"
    );
    assert_eq!(
        assessment.slope_semantic_role().as_str(),
        "directed-normal/minimum-slope"
    );
    assert_eq!(
        assessment.direction().components(),
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]
    );
    let recovered = assessment.recovered_gradient().components();
    for (actual, expected) in recovered.into_iter().zip([1.0, 2.0, 2.0]) {
        assert_close(actual, expected);
    }
    for residual in assessment.projection_residual().components() {
        assert_close(residual, 0.0);
    }
    assert_close(assessment.projection_residual_norm(), 0.0);
    assert_close(assessment.recovered_slope(), 3.0);
    assert_eq!(assessment.minimum_slope().value(), 2.5);
    assert_close(assessment.slope_slack(), 0.5);
    assert_close(assessment.slope_violation(), 0.0);
    assert_eq!(assessment.slope_active_state(), BoundActiveState::Inactive);
    assert_eq!(assessment.direction_loss(), None);
    assert_eq!(assessment.slope_loss(), None);
    assert_eq!(assessment.polarity_resolution_source_id(), None);
    assert_close(fit.report().field_energy().unwrap(), 0.0);
    assert_close(fit.report().total_objective().unwrap(), 0.0);
}

#[test]
fn axis_aligned_normal_problem_size_counts_the_retained_projection_block() {
    let mut problem = builder();
    add_affine_values(&mut problem, [0.0, 0.0, 1.0]);
    problem
        .add(
            DirectedNormalObservation::try_new(
                SourceId::new("axis-aligned-normal"),
                point(0.25, -0.5, 0.75),
                vector(0.0, 0.0, 1.0),
                MinimumNormalSlope::try_new(0.5).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();

    let fit = problem.build().unwrap().fit().unwrap();
    assert_eq!(fit.report().problem_size().scalar_hard_relations(), 8);
    assert_eq!(fit.report().directed_normals().len(), 1);
}

#[test]
fn soft_directed_normal_keeps_direction_and_slope_losses_independent() {
    let location = point(0.25, -0.5, 0.75);
    let mut problem = builder();
    add_affine_values(&mut problem, [1.0, 2.0, 0.0]);
    problem
        .add(GradientObservation::new(
            SourceId::new("fixed-gradient"),
            location,
            vector(1.0, 2.0, 0.0),
        ))
        .unwrap();
    problem
        .add(
            DirectedNormalObservation::try_with_enforcement(
                SourceId::new("soft-normal"),
                location,
                vector(0.0, 0.0, 4.0),
                NormalDirectionEnforcement::with_quadratic_penalty(
                    QuadraticPenalty::try_new(2.0).unwrap(),
                ),
                MinimumNormalSlope::try_new(0.5).unwrap(),
                MinimumNormalSlopeEnforcement::with_linear_violation_penalty(
                    LinearViolationPenalty::try_new(3.0).unwrap(),
                ),
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let fit = problem.build().unwrap().fit().unwrap();
    let assessment = &fit.report().directed_normals()[0];
    for (actual, expected) in assessment
        .projection_residual()
        .components()
        .into_iter()
        .zip([1.0, 2.0, 0.0])
    {
        assert_close(actual, expected);
    }
    assert_close(assessment.projection_residual_norm(), 5.0_f64.sqrt());
    assert_close(assessment.recovered_slope(), 0.0);
    assert_close(assessment.slope_slack(), 0.0);
    assert_close(assessment.slope_violation(), 0.5);
    assert_eq!(assessment.slope_active_state(), BoundActiveState::Active);
    assert_eq!(
        assessment.direction_quadratic_penalty().unwrap().weight(),
        2.0
    );
    assert_eq!(assessment.direction_standard_deviation(), None);
    assert_eq!(
        assessment
            .slope_linear_violation_penalty()
            .unwrap()
            .weight(),
        3.0
    );
    assert_close(assessment.direction_loss().unwrap(), 5.0);
    assert_close(assessment.slope_loss().unwrap(), 1.5);
    assert_close(fit.report().field_energy().unwrap(), 0.0);
    assert_close(fit.report().total_objective().unwrap(), 6.5);
}

#[test]
fn statistical_direction_and_quadratic_slope_channels_restore_their_own_losses() {
    let location = point(0.25, -0.5, 0.75);
    let mut problem = builder();
    add_affine_values(&mut problem, [1.0, 2.0, 0.0]);
    problem
        .add(GradientObservation::new(
            SourceId::new("statistical-fixed-gradient"),
            location,
            vector(1.0, 2.0, 0.0),
        ))
        .unwrap();
    problem
        .add(
            DirectedNormalObservation::try_with_enforcement(
                SourceId::new("statistical-normal"),
                location,
                vector(0.0, 0.0, 1.0),
                NormalDirectionEnforcement::with_standard_deviation(
                    StandardDeviation::try_new(0.5).unwrap(),
                ),
                MinimumNormalSlope::try_new(0.5).unwrap(),
                MinimumNormalSlopeEnforcement::with_quadratic_penalty(
                    QuadraticPenalty::try_new(4.0).unwrap(),
                ),
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let fit = problem.build().unwrap().fit().unwrap();
    let assessment = &fit.report().directed_normals()[0];
    assert_eq!(
        assessment.direction_standard_deviation().unwrap().value(),
        0.5
    );
    assert_eq!(assessment.slope_quadratic_penalty().unwrap().weight(), 4.0);
    assert_close(assessment.direction_loss().unwrap(), 10.0);
    assert_close(assessment.slope_loss().unwrap(), 0.5);
    assert_close(fit.report().total_objective().unwrap(), 10.5);
}

#[test]
fn hard_zero_gradient_does_not_satisfy_a_directed_normal() {
    let location = point(0.0, 0.0, 0.0);
    let mut problem = builder();
    add_affine_values(&mut problem, [0.0, 0.0, 0.0]);
    problem
        .add(GradientObservation::new(
            SourceId::new("zero-gradient"),
            location,
            vector(0.0, 0.0, 0.0),
        ))
        .unwrap();
    problem
        .add(
            DirectedNormalObservation::try_new(
                SourceId::new("requires-positive-slope"),
                location,
                vector(0.0, 0.0, 1.0),
                MinimumNormalSlope::try_new(0.25).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(!failure.report().direct_input_conflicts().is_empty());
    assert!(failure.report().attempts().is_empty());
}

#[test]
fn resolved_axial_normal_matches_directed_observables_and_retains_resolution_provenance() {
    let location = point(0.25, -0.5, 0.75);
    let slope = MinimumNormalSlope::try_new(2.5).unwrap();
    let mut directed_problem = builder();
    add_affine_values(&mut directed_problem, [1.0, 2.0, 2.0]);
    add_affine_shared_level(&mut directed_problem);
    directed_problem
        .add(
            DirectedNormalObservation::try_new(
                SourceId::new("normal"),
                location,
                vector(1.0, 2.0, 2.0),
                slope,
            )
            .unwrap(),
        )
        .unwrap();
    let directed = directed_problem.build().unwrap().fit().unwrap();

    let mut axial_problem = builder();
    add_affine_values(&mut axial_problem, [1.0, 2.0, 2.0]);
    add_affine_shared_level(&mut axial_problem);
    axial_problem
        .add(PolarityResolution::new(
            SourceId::new("field-resolution"),
            SourceId::new("normal"),
            PolaritySelection::AgainstInputAxis,
        ))
        .unwrap();
    axial_problem
        .add(
            AxialNormalObservation::try_new(
                SourceId::new("normal"),
                location,
                vector(-1.0, -2.0, -2.0),
                slope,
            )
            .unwrap(),
        )
        .unwrap();
    let axial = axial_problem.build().unwrap().fit().unwrap();

    let directed_assessment = &directed.report().directed_normals()[0];
    let axial_assessment = &axial.report().directed_normals()[0];
    assert_eq!(
        directed_assessment.direction(),
        axial_assessment.direction()
    );
    assert_close(
        directed_assessment.recovered_slope(),
        axial_assessment.recovered_slope(),
    );
    assert_eq!(
        axial_assessment.polarity_resolution_source_id(),
        Some(&SourceId::new("field-resolution"))
    );
    assert_eq!(
        axial_assessment.polarity_selection(),
        Some(PolaritySelection::AgainstInputAxis)
    );
    assert_close(
        directed
            .model()
            .shared_level_value(&GroupId::new("normal-level"))
            .unwrap(),
        axial
            .model()
            .shared_level_value(&GroupId::new("normal-level"))
            .unwrap(),
    );
    for query in [point(0.1, 0.2, 0.3), point(-0.5, 0.75, -0.25)] {
        let left = directed.model().evaluate(query).unwrap();
        let right = axial.model().evaluate(query).unwrap();
        assert_close(left.value(), right.value());
        for (left, right) in left
            .gradient()
            .components()
            .into_iter()
            .zip(right.gradient().components())
        {
            assert_close(left, right);
        }
    }
}

#[derive(Clone, Copy)]
enum FrameCase {
    Original,
    Rotated,
    Reflected,
}

fn orthogonal_transform(case: FrameCase, [x, y, z]: [f64; 3]) -> [f64; 3] {
    match case {
        FrameCase::Original => [x, y, z],
        FrameCase::Rotated => [-y, x, z],
        FrameCase::Reflected => [y, x, z],
    }
}

fn transformed_point(case: FrameCase, components: [f64; 3]) -> [f64; 3] {
    if matches!(case, FrameCase::Original) {
        components
    } else {
        let [x, y, z] = orthogonal_transform(case, components);
        [3.0 * x + 10.0, 3.0 * y - 4.0, 3.0 * z + 2.0]
    }
}

#[test]
fn normal_projection_uses_the_physical_euclidean_frame_not_kernel_anisotropy() {
    for case in [
        FrameCase::Original,
        FrameCase::Rotated,
        FrameCase::Reflected,
    ] {
        let (labels, handedness, scale) = match case {
            FrameCase::Original => (["x", "y", "z"], Handedness::Right, 1.0),
            FrameCase::Rotated => (["-y", "x", "z"], Handedness::Right, 3.0),
            FrameCase::Reflected => (["y", "x", "z"], Handedness::Left, 3.0),
        };
        let frame = InputCoordinateFrame::try_new(
            labels,
            handedness,
            LengthUnitLabel::new("transformed-length"),
        )
        .unwrap();
        let mut problem = ProblemBuilder::new(frame, FieldUnitLabel::new("field"));
        let original_gradient = [1.0, 2.0, 2.0];
        for (index, original) in [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, 1.0, 0.5],
        ]
        .into_iter()
        .enumerate()
        {
            let transformed = transformed_point(case, original);
            let value = 4.0
                + original_gradient[0] * original[0]
                + original_gradient[1] * original[1]
                + original_gradient[2] * original[2];
            problem
                .add(
                    georbf::observation::FieldValueObservation::try_new(
                        SourceId::new(format!("frame-value-{index}")),
                        point(transformed[0], transformed[1], transformed[2]),
                        value,
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        let support = transformed_point(case, [0.25, -0.5, 0.75]);
        let direction = orthogonal_transform(case, original_gradient);
        problem
            .add(
                DirectedNormalObservation::try_new(
                    SourceId::new("frame-normal"),
                    point(support[0], support[1], support[2]),
                    vector(direction[0], direction[1], direction[2]),
                    MinimumNormalSlope::try_new(2.5 / scale).unwrap(),
                )
                .unwrap(),
            )
            .unwrap();
        let metric = if matches!(case, FrameCase::Original) {
            [[4.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.25]]
        } else {
            [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 0.25]]
        };
        problem
            .set_global_anisotropy_metric(GlobalAnisotropyMetric::try_from_matrix(metric).unwrap())
            .unwrap();

        let fit = problem.build().unwrap().fit().unwrap();
        let normal = &fit.report().directed_normals()[0];
        assert_close(normal.projection_residual_norm(), 0.0);
        assert_close(normal.recovered_slope(), 3.0 / scale);
        let expected_gradient =
            orthogonal_transform(case, original_gradient).map(|component| component / scale);
        for (actual, expected) in normal
            .recovered_gradient()
            .components()
            .into_iter()
            .zip(expected_gradient)
        {
            assert_close(actual, expected);
        }
    }
}

#[test]
fn unresolved_semantics_has_priority_while_direct_conflict_evidence_is_retained() {
    let location = point(0.0, 0.0, 0.0);
    let mut problem = builder();
    add_affine_values(&mut problem, [0.0, 0.0, 0.0]);
    problem
        .add(GradientObservation::new(
            SourceId::new("priority-zero-gradient"),
            location,
            vector(0.0, 0.0, 0.0),
        ))
        .unwrap();
    problem
        .add(
            DirectedNormalObservation::try_new(
                SourceId::new("priority-directed"),
                location,
                vector(0.0, 0.0, 1.0),
                MinimumNormalSlope::try_new(0.25).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(
            AxialNormalObservation::try_new(
                SourceId::new("priority-unresolved"),
                location,
                vector(1.0, 0.0, 0.0),
                MinimumNormalSlope::try_new(0.25).unwrap(),
            )
            .unwrap(),
        )
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::UnresolvedSemantics);
    assert_eq!(failure.report().unresolved_axial_normals().len(), 1);
    assert!(!failure.report().direct_input_conflicts().is_empty());
    assert!(failure.report().attempts().is_empty());
}

#[test]
fn nonlocal_normal_geometry_requires_a_validated_infeasibility_certificate() {
    let location = point(0.0, 0.0, 0.0);
    let mut problem = builder();
    add_affine_values(&mut problem, [1.0, 0.0, 0.0]);
    for (source, direction) in [
        ("normal-x", [1.0, 0.0, 0.0]),
        ("normal-y", [0.0, 1.0, 0.0]),
        ("normal-negative-xy", [-1.0, -1.0, 0.0]),
    ] {
        problem
            .add(
                DirectedNormalObservation::try_with_enforcement(
                    SourceId::new(source),
                    location,
                    vector(direction[0], direction[1], direction[2]),
                    NormalDirectionEnforcement::with_quadratic_penalty(
                        QuadraticPenalty::try_new(1.0).unwrap(),
                    ),
                    MinimumNormalSlope::try_new(1.0).unwrap(),
                    MinimumNormalSlopeEnforcement::hard(),
                )
                .unwrap(),
            )
            .unwrap();
    }
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    eprintln!("{:#?}", failure.report());
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::InfeasibleProblem);
    assert!(failure.report().direct_input_conflicts().is_empty());
    let certificate = failure
        .report()
        .infeasibility_certificate()
        .expect("general affine infeasibility requires an independently checked ray");
    assert!(certificate.finite());
    assert!(certificate.backend_invoked());
}
