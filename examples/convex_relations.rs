use std::error::Error;

use georbf::diagnostics::SolveAttemptTermination;
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::{FieldEnergyNormalization, KernelKind};
use georbf::observation::{
    AxialNormalObservation, CovarianceGroupBuilder, CovarianceMatrix, DirectedNormalObservation,
    FieldValueObservation, GradientObservation, MinimumNormalSlope, MinimumNormalSlopeEnforcement,
    NormalDirectionEnforcement, QuadraticPenalty, TangentDirectionObservation,
};
use georbf::relation::HorizonBuilder;
use georbf::relation::{
    DirectionalDerivativeInterval, FieldLevelOrder, FieldSeparationInterval,
    FieldSeparationViolationPenalty, FieldValueBound, LinearViolationPenalty, MinimumFieldOffset,
    MinimumFieldSeparation, PointToLevelSetRelation, PointToLevelSetSide, PolarityResolution,
    PolaritySelection, StratigraphicFieldDirection, YoungerThan,
};
use georbf::{GroupId, Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("manufactured example points are finite")
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).expect("manufactured example vectors are finite")
}

fn planar_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    x + 2.0 * y + 3.0 * z
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-8 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

/// Runs the complete public v0.2.0 Convex Relations preview workflow.
pub fn run() -> Result<(), Box<dyn Error>> {
    // Manufactured truth: f(x, y, z) = x + 2y + 3z.
    let frame = InputCoordinateFrame::try_new(
        ["east", "north", "elevation"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )?;
    let mut problem = ProblemBuilder::new(frame, FieldUnitLabel::new("stratigraphic-unit"));
    problem.set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)?;

    for (source, location) in [
        ("hard-value/origin", point(0.0, 0.0, 0.0)),
        ("hard-value/east", point(1.0, 0.0, 0.0)),
        ("hard-value/north", point(0.0, 1.0, 0.0)),
        ("hard-value/up", point(0.0, 0.0, 1.0)),
        ("hard-value/diagonal", point(1.0, 1.0, 0.0)),
    ] {
        problem.add(FieldValueObservation::try_new(
            SourceId::new(source),
            location,
            planar_value(location),
        )?)?;
    }

    let quadratic = QuadraticPenalty::try_new(2.0)?;
    let linear = LinearViolationPenalty::try_new(2.0)?;
    problem.add(FieldValueObservation::try_with_quadratic_penalty(
        SourceId::new("soft-value"),
        point(0.5, 0.5, 0.5),
        3.0,
        quadratic,
    )?)?;
    problem.add(GradientObservation::new(
        SourceId::new("hard-gradient"),
        point(0.25, -0.5, 0.75),
        vector(1.0, 2.0, 3.0),
    ))?;
    problem.add(GradientObservation::with_quadratic_penalty(
        SourceId::new("soft-gradient"),
        point(-0.5, 0.25, 0.5),
        vector(1.0, 2.0, 3.0),
        quadratic,
    ))?;
    problem.add(TangentDirectionObservation::try_new(
        SourceId::new("hard-tangent"),
        point(-0.75, 0.25, 0.5),
        vector(2.0, -1.0, 0.0),
    )?)?;
    problem.add(TangentDirectionObservation::try_with_quadratic_penalty(
        SourceId::new("soft-tangent"),
        point(0.75, -0.25, 0.5),
        vector(3.0, 0.0, -1.0),
        quadratic,
    )?)?;

    let reference_id = GroupId::new("older-horizon");
    let mut reference = HorizonBuilder::new(reference_id.clone());
    reference.add_member(SourceId::new("older/origin"), point(0.0, 0.0, 0.0))?;
    reference.add_member(SourceId::new("older/second"), point(2.0, -1.0, 0.0))?;
    problem.add(reference.build()?)?;

    let target_id = GroupId::new("younger-horizon");
    let mut target = HorizonBuilder::new(target_id.clone());
    target.add_member(SourceId::new("younger/diagonal"), point(1.0, 1.0, 0.0))?;
    target.add_member(SourceId::new("younger/up"), point(0.0, 0.0, 1.0))?;
    problem.add(target.build()?)?;

    problem.add(FieldValueBound::try_interval(
        SourceId::new("hard-field-bound"),
        point(0.0, 2.0, 0.0),
        3.5,
        4.5,
    )?)?;
    problem.add(FieldValueBound::try_interval_with_violation_penalties(
        SourceId::new("soft-field-bound"),
        point(0.0, 0.0, 2.0),
        5.5,
        quadratic.into(),
        6.5,
        linear.into(),
    )?)?;
    problem.add(DirectionalDerivativeInterval::try_interval(
        SourceId::new("hard-derivative-interval"),
        point(0.25, -0.5, 0.75),
        vector(1.0, 0.0, 0.0),
        0.9,
        1.1,
    )?)?;
    problem.add(
        DirectionalDerivativeInterval::try_interval_with_violation_penalties(
            SourceId::new("soft-derivative-interval"),
            point(-0.25, 0.5, 0.25),
            vector(0.0, 1.0, 0.0),
            1.9,
            quadratic.into(),
            2.1,
            linear.into(),
        )?,
    )?;

    let minimum_separation = MinimumFieldSeparation::try_new(2.0)?;
    problem.add(YoungerThan::hard(
        SourceId::new("hard-horizon-order"),
        target_id.clone(),
        reference_id.clone(),
        minimum_separation,
    ))?;
    problem.add(YoungerThan::with_quadratic_penalty(
        SourceId::new("soft-horizon-order"),
        target_id.clone(),
        reference_id.clone(),
        minimum_separation,
        quadratic,
    ))?;
    problem.add(FieldLevelOrder::hard(
        SourceId::new("hard-field-level-order"),
        reference_id.clone(),
        target_id.clone(),
    ))?;
    problem.add(FieldSeparationInterval::try_hard(
        SourceId::new("hard-field-separation"),
        reference_id.clone(),
        target_id.clone(),
        2.5,
        3.5,
    )?)?;
    problem.add(FieldSeparationInterval::try_with_violation_penalties(
        SourceId::new("soft-field-separation"),
        reference_id.clone(),
        target_id.clone(),
        2.75,
        FieldSeparationViolationPenalty::Quadratic(quadratic),
        3.25,
        FieldSeparationViolationPenalty::Linear(linear),
    )?)?;
    problem.add(PointToLevelSetRelation::hard(
        SourceId::new("hard-point-side"),
        point(1.0, 0.0, 1.0),
        target_id.clone(),
        PointToLevelSetSide::Increasing,
        MinimumFieldOffset::try_new(0.5)?,
    ))?;
    problem.add(PointToLevelSetRelation::with_quadratic_penalty(
        SourceId::new("soft-point-side"),
        point(0.0, 1.0, 1.0),
        target_id.clone(),
        PointToLevelSetSide::Increasing,
        MinimumFieldOffset::try_new(1.5)?,
        quadratic,
    ))?;

    let normal = vector(1.0, 2.0, 3.0);
    let minimum_slope = MinimumNormalSlope::try_new(3.0)?;
    problem.add(DirectedNormalObservation::try_new(
        SourceId::new("hard-directed-normal"),
        point(0.2, -0.3, 0.4),
        normal,
        minimum_slope,
    )?)?;
    problem.add(PolarityResolution::new(
        SourceId::new("soft-axial-resolution"),
        SourceId::new("soft-axial-normal"),
        PolaritySelection::AlongInputAxis,
    ))?;
    problem.add(AxialNormalObservation::try_with_enforcement(
        SourceId::new("soft-axial-normal"),
        point(-0.4, 0.2, 0.1),
        normal,
        NormalDirectionEnforcement::with_quadratic_penalty(quadratic),
        minimum_slope,
        MinimumNormalSlopeEnforcement::with_linear_violation_penalty(linear),
    )?)?;

    let mut covariance = CovarianceGroupBuilder::new(GroupId::new("field-covariance"));
    covariance.add_field_value_member(
        SourceId::new("covariance/east"),
        point(1.0, 0.0, 0.0),
        1.0,
    )?;
    covariance.add_field_value_member(
        SourceId::new("covariance/north"),
        point(0.0, 1.0, 0.0),
        2.0,
    )?;
    problem.add(covariance.build(CovarianceMatrix::try_new([[1.0, 0.25], [0.25, 2.0]])?)?)?;
    problem.set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0)?)?;

    let snapshot = problem.build()?;
    assert_eq!(snapshot.resolved_kernel().kind(), KernelKind::Cubic);
    let success = snapshot.fit()?;
    let report = success.report();
    let acceptance = report
        .canonical_acceptance()
        .expect("a successful fit records canonical acceptance");
    assert!(acceptance.accepted() && acceptance.provenance_verified());
    assert!(acceptance.objective_verified());
    assert!(
        acceptance
            .hard_affine_inequality_violation_max()
            .unwrap_or(0.0)
            <= 1.0e-8
    );
    assert!(acceptance.scaling_round_trip_error().unwrap_or(0.0) <= 1.0e-11);
    assert!(acceptance.reduction_round_trip_error().unwrap_or(0.0) <= 1.0e-11);
    assert_eq!(report.field_value_bounds().len(), 4);
    assert_eq!(report.directional_derivative_intervals().len(), 4);
    assert_eq!(report.shared_level_set_relations().len(), 3);
    assert_eq!(report.field_separation_intervals().len(), 4);
    assert_eq!(report.point_to_level_set_relations().len(), 2);
    assert_eq!(report.directed_normals().len(), 2);
    assert_eq!(report.soft_field_values().len(), 1);
    assert_eq!(report.soft_gradients().len(), 1);
    assert_eq!(report.soft_tangents().len(), 1);
    assert_eq!(report.covariance_groups().len(), 1);
    assert_eq!(report.shared_level_values().len(), 2);
    assert!(matches!(
        report
            .attempts()
            .last()
            .map(|attempt| attempt.termination()),
        Some(
            SolveAttemptTermination::CandidateProduced
                | SolveAttemptTermination::ReducedAccuracyCandidateProduced
        )
    ));

    for group_id in [&reference_id, &target_id] {
        assert!(success.model().shared_level_value(group_id).is_some());
    }
    let queries = [
        point(0.75, -1.0, 2.0),
        point(-0.5, 0.25, -0.75),
        point(1.5, 2.0, 0.5),
    ];
    let batch = success.model().evaluate_batch(&queries)?;
    for ((location, sample), single) in queries
        .into_iter()
        .zip(batch)
        .zip(queries.map(|location| success.model().evaluate(location)))
    {
        assert_eq!(sample, single?);
        assert_close(sample.value(), planar_value(location));
        for (actual, expected) in sample
            .gradient()
            .components()
            .into_iter()
            .zip([1.0, 2.0, 3.0])
        {
            assert_close(actual, expected);
        }
    }

    Ok(())
}

/// Runs the release quantity smoke without establishing a timing SLA.
pub fn run_smoke() -> Result<(), Box<dyn Error>> {
    let frame = InputCoordinateFrame::try_new(
        ["east", "north", "elevation"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )?;
    let mut problem = ProblemBuilder::new(frame, FieldUnitLabel::new("smoke-field"));
    for (index, location) in [
        point(-1.0, -1.0, -1.0),
        point(1.0, -1.0, -1.0),
        point(-1.0, 1.0, -1.0),
        point(-1.0, -1.0, 1.0),
        point(1.0, 1.0, 0.5),
    ]
    .into_iter()
    .enumerate()
    {
        problem.add(FieldValueObservation::try_new(
            SourceId::new(format!("smoke/value-{index:03}")),
            location,
            planar_value(location),
        )?)?;
    }
    for index in 0..507 {
        let location = point(
            f64::from(index % 13) - 6.0,
            f64::from((index / 13) % 13) - 6.0,
            f64::from(index / 169) - 1.0,
        );
        problem.add(FieldValueBound::try_upper(
            SourceId::new(format!("smoke/bound-{index:03}")),
            location,
            planar_value(location) + 1.0,
        )?)?;
    }

    let success = problem.build()?.fit()?;
    assert_eq!(success.report().problem_size().scalar_hard_relations(), 512);
    assert_eq!(success.report().field_value_bounds().len(), 507);
    let queries = (0..10_000)
        .map(|index| {
            let x = f64::from(index % 101) / 10.0 - 5.0;
            let y = f64::from((index / 101) % 97) / 10.0 - 4.0;
            let z = f64::from(index % 29) / 20.0 - 0.5;
            point(x, y, z)
        })
        .collect::<Vec<_>>();
    let samples = success.model().evaluate_batch(&queries)?;
    assert_eq!(samples.len(), 10_000);
    for index in [0, 4_999, 9_999] {
        assert_eq!(samples[index], success.model().evaluate(queries[index])?);
    }

    Ok(())
}

#[cfg(not(test))]
fn main() -> Result<(), Box<dyn Error>> {
    run()
}
