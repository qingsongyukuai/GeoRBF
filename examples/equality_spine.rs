use std::collections::BTreeSet;
use std::error::Error;

use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::KernelKind;
use georbf::observation::{
    FieldValueObservation, GradientObservation, TangentDirectionObservation,
};
use georbf::relation::{AdditiveFieldGauge, HorizonBuilder};
use georbf::{GroupId, Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the manufactured point is finite")
}

fn planar_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    3.0 + 0.5 * x - 0.25 * y + z
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-9 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

/// Runs the complete public v0.1.0 Equality Spine workflow.
pub fn run() -> Result<(), Box<dyn Error>> {
    let frame = InputCoordinateFrame::try_new(
        ["east", "north", "elevation"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )?;
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("stratigraphic-unit"));

    let horizon_id = GroupId::new("planar-horizon");
    let mut horizon = HorizonBuilder::new(horizon_id.clone());
    for (source, location) in [
        ("horizon/left", point(-2.0, 0.0, 1.0)),
        ("horizon/right", point(2.0, 0.0, -1.0)),
        ("horizon/north", point(0.0, 4.0, 1.0)),
    ] {
        horizon.add_member(SourceId::new(source), location)?;
    }
    builder.add(horizon.build()?)?;
    builder.add(AdditiveFieldGauge::at_level_set(
        SourceId::new("gauge/planar-horizon"),
        horizon_id.clone(),
        3.0,
    )?)?;

    let derivative_location = point(0.5, -0.5, 0.25);
    builder.add(GradientObservation::new(
        SourceId::new("gradient/complete"),
        derivative_location,
        Vector3::try_new(0.5, -0.25, 1.0)?,
    ))?;
    builder.add(TangentDirectionObservation::try_new(
        SourceId::new("tangent/strike"),
        point(-0.75, 0.25, 0.5),
        Vector3::try_new(2.0, 0.0, -1.0)?,
    )?)?;

    let absolute_location = point(0.0, 0.0, 2.0);
    builder.add(FieldValueObservation::try_new(
        SourceId::new("field-value/absolute"),
        absolute_location,
        planar_value(absolute_location),
    )?)?;

    let snapshot = builder.build()?;
    assert_eq!(snapshot.horizon_count(), 1);
    assert_eq!(snapshot.resolved_kernel().kind(), KernelKind::Cubic);

    let success = snapshot.fit()?;
    let report = success.report();
    assert_eq!(report.resolved_kernel().kind(), KernelKind::Cubic);
    assert_eq!(report.numerical_policy().as_str(), "georbf-v1");
    assert_eq!(
        report
            .backend_fingerprint()
            .expect("a successful fit records its backend")
            .crate_name(),
        "faer"
    );
    assert!(
        report
            .canonical_acceptance()
            .is_some_and(|acceptance| acceptance.accepted() && acceptance.provenance_verified())
    );
    assert_close(
        success
            .model()
            .shared_level_value(&horizon_id)
            .expect("the recovered horizon latent is public by GroupId"),
        3.0,
    );

    let relation_sources = report
        .hard_relations()
        .iter()
        .map(|relation| relation.source_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        relation_sources,
        BTreeSet::from([
            "field-value/absolute",
            "gauge/planar-horizon",
            "gradient/complete",
            "horizon/left",
            "horizon/north",
            "horizon/right",
            "tangent/strike",
        ])
    );
    assert_eq!(report.hard_relations().len(), 9);
    assert!(report.hard_relations().iter().all(|relation| {
        !relation.semantic_role().as_str().is_empty()
            && relation.residual().abs() <= relation.tolerance()
    }));

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
            .zip([0.5, -0.25, 1.0])
        {
            assert_close(actual, expected);
        }
    }

    Ok(())
}

#[cfg(not(test))]
fn main() -> Result<(), Box<dyn Error>> {
    run()
}
