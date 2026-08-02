use georbf::diagnostics::ResidualDimension;
use georbf::geometry::{
    FieldUnitLabel, GeometryError, Handedness, InputCoordinateFrame, LengthUnitLabel, Point3,
    Vector3,
};
use georbf::observation::{FieldValueObservation, ObservationError, TangentDirectionObservation};
use georbf::problem::AddError;
use georbf::{ProblemBuilder, SourceId};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the manufactured point is finite")
}

fn vector(x: f64, y: f64, z: f64) -> Vector3 {
    Vector3::try_new(x, y, z).expect("the manufactured vector is finite")
}

fn problem_builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .expect("the manufactured frame is valid"),
        FieldUnitLabel::new("field-unit"),
    )
}

fn affine_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    4.0 + 2.0 * x - y + 0.5 * z
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-9 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn fit_affine_with_tangent(
    tangent_direction: [f64; 3],
    reverse_input_order: bool,
) -> georbf::fit::FitSuccess {
    let tangent = TangentDirectionObservation::try_new(
        SourceId::new("equivalent-tangent"),
        point(0.25, -0.5, 0.75),
        vector(
            tangent_direction[0],
            tangent_direction[1],
            tangent_direction[2],
        ),
    )
    .unwrap();
    let mut supports = [
        point(-1.0, -1.0, -1.0),
        point(1.0, -1.0, -1.0),
        point(-1.0, 1.0, -1.0),
        point(-1.0, -1.0, 1.0),
        point(1.0, 1.0, 0.5),
    ];
    if reverse_input_order {
        supports.reverse();
    }
    let mut builder = problem_builder();
    if reverse_input_order {
        builder.add(tangent.clone()).unwrap();
    }
    for location in supports {
        let [x, y, z] = location.components();
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("value/{x}/{y}/{z}")),
                    location,
                    affine_value(location),
                )
                .unwrap(),
            )
            .unwrap();
    }
    if !reverse_input_order {
        builder.add(tangent).unwrap();
    }
    builder.build().unwrap().fit().unwrap()
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

fn transform_point(case: FrameCase, components: [f64; 3]) -> [f64; 3] {
    if matches!(case, FrameCase::Original) {
        return components;
    }
    let [x, y, z] = orthogonal_transform(case, components);
    [3.0 * x + 10.0, 3.0 * y - 4.0, 3.0 * z + 2.0]
}

fn transform_gradient(case: FrameCase, components: [f64; 3]) -> [f64; 3] {
    if matches!(case, FrameCase::Original) {
        components
    } else {
        orthogonal_transform(case, components).map(|component| component / 3.0)
    }
}

fn fit_affine_in_frame(case: FrameCase) -> georbf::fit::FitSuccess {
    let (axis_labels, handedness, length_unit) = match case {
        FrameCase::Original => (["east", "north", "elevation"], Handedness::Right, "m"),
        FrameCase::Rotated => (
            ["rotated-north", "rotated-east", "rotated-elevation"],
            Handedness::Right,
            "scaled-m",
        ),
        FrameCase::Reflected => (
            ["reflected-north", "reflected-east", "reflected-elevation"],
            Handedness::Left,
            "scaled-m",
        ),
    };
    let frame =
        InputCoordinateFrame::try_new(axis_labels, handedness, LengthUnitLabel::new(length_unit))
            .unwrap();
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("field-unit"));
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
        let original = point(support[0], support[1], support[2]);
        let transformed = transform_point(case, support);
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("frame-value-{index}")),
                    point(transformed[0], transformed[1], transformed[2]),
                    affine_value(original),
                )
                .unwrap(),
            )
            .unwrap();
    }
    let tangent_support = [0.25, -0.5, 0.75];
    let tangent_direction = orthogonal_transform(case, [1.0, 2.0, 0.0]);
    let tangent_support = transform_point(case, tangent_support);
    builder
        .add(
            TangentDirectionObservation::try_new(
                SourceId::new("frame-tangent"),
                point(tangent_support[0], tangent_support[1], tangent_support[2]),
                vector(
                    tangent_direction[0],
                    tangent_direction[1],
                    tangent_direction[2],
                ),
            )
            .unwrap(),
        )
        .unwrap();
    builder.build().unwrap().fit().unwrap()
}

#[test]
fn checked_tangent_directions_are_unit_axial_observations() {
    let location = point(1.0, -2.0, 0.5);
    let tangent = TangentDirectionObservation::try_new(
        SourceId::new("tangent"),
        location,
        vector(-6.0, 8.0, -0.0),
    )
    .expect("a finite nonzero direction is valid");
    let reversed = TangentDirectionObservation::try_new(
        SourceId::new("reversed"),
        location,
        vector(6.0, -8.0, 0.0),
    )
    .expect("the opposite direction denotes the same tangent axis");

    assert_eq!(tangent.location(), location);
    assert_eq!(tangent.direction().components(), [0.6, -0.8, 0.0]);
    assert_eq!(reversed.direction(), tangent.direction());
    assert_eq!(tangent.source_id().as_str(), "tangent");
    assert!(
        tangent
            .direction()
            .components()
            .iter()
            .filter(|component| **component == 0.0)
            .all(|component| !component.is_sign_negative())
    );

    assert!(matches!(
        TangentDirectionObservation::try_new(
            SourceId::new("zero"),
            location,
            vector(0.0, -0.0, 0.0),
        ),
        Err(ObservationError::ZeroTangentDirection)
    ));
    assert!(matches!(
        Vector3::try_new(1.0, f64::INFINITY, 0.0),
        Err(GeometryError::NonFiniteVector { axis: 1 })
    ));

    for direction in [
        vector(f64::MAX, -f64::MAX, 0.0),
        vector(f64::MIN_POSITIVE, f64::MIN_POSITIVE, 0.0),
        vector(f64::from_bits(1), 0.0, 0.0),
    ] {
        let unit =
            TangentDirectionObservation::try_new(SourceId::new("extreme"), location, direction)
                .expect("every finite nonzero magnitude normalizes without overflow or underflow")
                .direction()
                .components();
        assert!(unit.iter().all(|component| component.is_finite()));
        assert_close(
            unit.into_iter()
                .map(|component| component * component)
                .sum::<f64>()
                .sqrt(),
            1.0,
        );
    }
}

#[test]
fn tangent_directions_enter_snapshots_with_atomic_source_identity() {
    let location = point(0.0, 0.0, 0.0);
    let mut builder = problem_builder();
    builder
        .add(
            TangentDirectionObservation::try_new(
                SourceId::new("tangent-a"),
                location,
                vector(1.0, 2.0, 0.0),
            )
            .unwrap(),
        )
        .expect("the first SourceId is unique");
    assert!(matches!(
        builder.add(
            TangentDirectionObservation::try_new(
                SourceId::new("tangent-a"),
                location,
                vector(0.0, 0.0, 1.0),
            )
            .unwrap(),
        ),
        Err(AddError::DuplicateSourceId { .. })
    ));
    builder
        .add(
            TangentDirectionObservation::try_new(
                SourceId::new("tangent-b"),
                location,
                vector(0.0, 1.0, 0.0),
            )
            .unwrap(),
        )
        .expect("a rejected duplicate did not mutate the builder");

    let snapshot = builder.build().expect("two hard tangents form a snapshot");
    assert_eq!(snapshot.observation_count(), 2);
    assert_eq!(snapshot.source_count(), 2);
}

#[test]
fn tangent_direction_fits_through_cubic_equality_and_is_physically_reverified() {
    let mut builder = problem_builder();
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
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("value-{index}")),
                    location,
                    affine_value(location),
                )
                .unwrap(),
            )
            .unwrap();
    }
    let tangent_location = point(0.25, -0.5, 0.75);
    builder
        .add(
            TangentDirectionObservation::try_new(
                SourceId::new("tangent"),
                tangent_location,
                vector(1.0, 2.0, 0.0),
            )
            .unwrap(),
        )
        .unwrap();

    let success = builder
        .build()
        .unwrap()
        .fit()
        .expect("the manufactured tangent relation is feasible");
    let query = point(0.2, -0.3, 0.4);
    let sample = success.model().evaluate(query).unwrap();
    assert_close(sample.value(), affine_value(query));
    for (actual, expected) in sample
        .gradient()
        .components()
        .into_iter()
        .zip([2.0, -1.0, 0.5])
    {
        assert_close(actual, expected);
    }

    let tangent_sample = success.model().evaluate(tangent_location).unwrap();
    let [gx, gy, _] = tangent_sample.gradient().components();
    let independently_recomputed = gx / 5.0_f64.sqrt() + 2.0 * gy / 5.0_f64.sqrt();
    assert_close(independently_recomputed, 0.0);

    let report = success.report();
    assert_eq!(report.problem_size().input_observations(), 6);
    assert_eq!(report.problem_size().scalar_hard_relations(), 6);
    let tangent_relation = report
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id().as_str() == "tangent")
        .expect("the tangent relation retains its caller source");
    assert_eq!(
        tangent_relation.semantic_role().as_str(),
        "tangent-direction-observation/directional-derivative"
    );
    assert_eq!(
        tangent_relation.dimension(),
        ResidualDimension::FieldValuePerLength
    );
    assert_eq!(tangent_relation.target(), 0.0);
    assert_close(tangent_relation.recovered_value(), independently_recomputed);
    assert!(tangent_relation.residual().abs() <= tangent_relation.tolerance());
    assert_close(report.field_energy().unwrap(), 0.0);
    assert!(
        report
            .canonical_acceptance()
            .expect("successful recovery retains canonical evidence")
            .provenance_verified()
    );
}

#[test]
fn zero_gradient_satisfies_tangent_directions_without_a_slope_assumption() {
    let origin = point(0.0, 0.0, 0.0);
    let mut builder = problem_builder();
    for (index, location, value) in [
        (0, origin, 0.0),
        (1, point(1.0, 0.0, 0.0), 1.0),
        (2, point(0.0, 1.0, 0.0), 2.0),
        (3, point(0.0, 0.0, 1.0), 3.0),
        (4, point(-1.0, -1.0, -1.0), 6.0),
    ] {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("stationary-value-{index}")),
                    location,
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    for (source, location, direction) in [
        ("tangent-x", origin, vector(1.0, 0.0, 0.0)),
        ("tangent-y", origin, vector(0.0, 1.0, 0.0)),
        ("tangent-z", origin, vector(0.0, 0.0, 1.0)),
    ] {
        builder
            .add(
                TangentDirectionObservation::try_new(SourceId::new(source), location, direction)
                    .unwrap(),
            )
            .unwrap();
    }

    let success = builder
        .build()
        .unwrap()
        .fit()
        .expect("zero gradient at the tangent support is feasible");
    let stationary = success.model().evaluate(origin).unwrap();
    assert_close(stationary.value(), 0.0);
    for component in stationary.gradient().components() {
        assert_close(component, 0.0);
    }
    assert_eq!(
        success
            .report()
            .hard_relations()
            .iter()
            .filter(|relation| relation.source_id().as_str().starts_with("tangent-"))
            .count(),
        3
    );
    assert!(
        success
            .report()
            .hard_relations()
            .iter()
            .all(|relation| relation.residual().abs() <= relation.tolerance())
    );
}

#[test]
fn opposite_scaled_tangents_and_input_permutations_have_one_canonical_semantics() {
    let baseline = fit_affine_with_tangent([1.0, 2.0, 0.0], false);
    let equivalent = fit_affine_with_tangent([-8.0, -16.0, -0.0], true);
    let query = point(0.75, -0.25, 1.5);
    let baseline_sample = baseline.model().evaluate(query).unwrap();
    let equivalent_sample = equivalent.model().evaluate(query).unwrap();

    assert_eq!(equivalent_sample, baseline_sample);
    assert_eq!(
        equivalent.report().problem_size(),
        baseline.report().problem_size()
    );
    assert_eq!(
        equivalent.report().hard_relations(),
        baseline.report().hard_relations()
    );
    assert_eq!(
        equivalent.report().field_energy(),
        baseline.report().field_energy()
    );
    assert_eq!(
        equivalent.report().total_objective(),
        baseline.report().total_objective()
    );
    assert_eq!(
        equivalent
            .report()
            .canonical_acceptance()
            .unwrap()
            .provenance_verified(),
        baseline
            .report()
            .canonical_acceptance()
            .unwrap()
            .provenance_verified()
    );
}

#[test]
fn tangent_semantics_are_covariant_under_rotation_reflection_and_uniform_scale() {
    let original = fit_affine_in_frame(FrameCase::Original);
    let query = [0.2, -0.3, 0.4];
    let original_sample = original
        .model()
        .evaluate(point(query[0], query[1], query[2]))
        .unwrap();
    let original_tangent = original
        .report()
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id().as_str() == "frame-tangent")
        .unwrap();

    for case in [FrameCase::Rotated, FrameCase::Reflected] {
        let transformed = fit_affine_in_frame(case);
        let transformed_query = transform_point(case, query);
        let transformed_sample = transformed
            .model()
            .evaluate(point(
                transformed_query[0],
                transformed_query[1],
                transformed_query[2],
            ))
            .unwrap();
        assert_close(transformed_sample.value(), original_sample.value());
        for (actual, expected) in
            transformed_sample
                .gradient()
                .components()
                .into_iter()
                .zip(transform_gradient(
                    case,
                    original_sample.gradient().components(),
                ))
        {
            assert_close(actual, expected);
        }

        let transformed_tangent = transformed
            .report()
            .hard_relations()
            .iter()
            .find(|relation| relation.source_id().as_str() == "frame-tangent")
            .unwrap();
        assert_eq!(
            transformed_tangent.semantic_role(),
            original_tangent.semantic_role()
        );
        assert_eq!(
            transformed_tangent.dimension(),
            original_tangent.dimension()
        );
        assert_close(
            3.0 * transformed_tangent.recovered_value(),
            original_tangent.recovered_value(),
        );
        assert_close(
            3.0 * transformed_tangent.characteristic_scale(),
            original_tangent.characteristic_scale(),
        );
        assert_close(
            3.0 * transformed_tangent.tolerance(),
            original_tangent.tolerance(),
        );
        assert!(transformed_tangent.residual().abs() <= transformed_tangent.tolerance());
        assert_eq!(
            transformed.report().problem_size(),
            original.report().problem_size()
        );
        assert!(
            transformed
                .report()
                .canonical_acceptance()
                .unwrap()
                .accepted()
        );
    }
}
