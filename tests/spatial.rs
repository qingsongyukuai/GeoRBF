use georbf::{
    average_nearest_neighbour_distance, bounds, closest_to_distance_index, constraints_to_points,
    distance_between_points, extremal_point_indices, farthest_from_other_set_index,
    farthest_neighbour_index, farthest_pair_indices, largest_distance_between_points,
    maximal_axial_variability_order, nearest_neighbour_index, nearest_neighbour_indices,
    spatial_metrics, Axis, Constraints, Inequality, Interface, Planar, Point, Polarity,
    SpatialError, Tangent,
};

const EPSILON: f64 = 2.0e-14;

fn point(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z).unwrap()
}

fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= EPSILON * scale,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn conversion_order_and_four_dimensional_distance_match_surfe() {
    let constraints = Constraints {
        inequalities: vec![Inequality::with_c(1.0, 0.0, 0.0, 10.0, 11.0).unwrap()],
        interfaces: vec![Interface::with_c(2.0, 0.0, 0.0, 20.0, 22.0).unwrap()],
        planars: vec![Planar::from_normal_with_c(3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 33.0).unwrap()],
        tangents: vec![Tangent::with_c(4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 44.0).unwrap()],
    };
    let points = constraints_to_points(&constraints);
    assert_eq!(
        points
            .iter()
            .map(|value| (value.x(), value.c()))
            .collect::<Vec<_>>(),
        [(1.0, 11.0), (2.0, 22.0), (3.0, 33.0), (4.0, 44.0)]
    );

    let first = Point::with_c(0.0, 0.0, 0.0, 0.0).unwrap();
    let second = Point::with_c(1.0, 2.0, 2.0, 2.0).unwrap();
    assert_close(distance_between_points(&first, &second), 13.0_f64.sqrt());
    assert_close(distance_between_points(&second, &first), 13.0_f64.sqrt());
    assert_eq!(distance_between_points(&first, &first).to_bits(), 0);
}

#[test]
fn nearest_neighbours_skip_zero_distance_and_break_ties_deterministically() {
    let query = point(0.0, 0.0, 0.0);
    let candidates = vec![
        point(1.0, 0.0, 0.0),
        point(-1.0, 0.0, 0.0),
        point(0.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0),
    ];

    assert_eq!(nearest_neighbour_index(&query, &candidates), Ok(0));
    assert_eq!(nearest_neighbour_indices(2, &query, &candidates), [0, 1]);
    assert_eq!(nearest_neighbour_indices(4, &query, &candidates), [0, 1, 3]);
    assert!(nearest_neighbour_indices(0, &query, &candidates).is_empty());
    assert!(nearest_neighbour_indices(3, &query, &[]).is_empty());

    assert_eq!(
        nearest_neighbour_index(&query, &[]),
        Err(SpatialError::EmptyPointSet)
    );
    assert_eq!(
        nearest_neighbour_index(&query, &[point(0.0, 0.0, 0.0)]),
        Err(SpatialError::NoNonzeroNeighbour)
    );
}

#[test]
fn farthest_and_target_distance_helpers_preserve_first_tie_indices() {
    let query = point(0.0, 0.0, 0.0);
    let candidates = vec![
        point(2.0, 0.0, 0.0),
        point(-2.0, 0.0, 0.0),
        point(1.0, 0.0, 0.0),
    ];
    assert_eq!(farthest_neighbour_index(&query, &candidates), Ok(0));
    assert_eq!(closest_to_distance_index(&query, &candidates, 1.5), Ok(0));

    let sources = vec![point(1.0, 0.0, 0.0), point(-1.0, 0.0, 0.0)];
    let references = vec![point(0.0, 0.0, 0.0)];
    assert_eq!(farthest_from_other_set_index(&sources, &references), Ok(0));

    assert_eq!(
        farthest_neighbour_index(&query, &[]),
        Err(SpatialError::EmptyPointSet)
    );
    assert_eq!(
        farthest_from_other_set_index(&[], &references),
        Err(SpatialError::EmptyPointSet)
    );
    assert_eq!(
        farthest_from_other_set_index(&sources, &[]),
        Err(SpatialError::EmptyPointSet)
    );
    assert_eq!(
        closest_to_distance_index(&query, &[], 1.0),
        Err(SpatialError::EmptyPointSet)
    );
    assert_eq!(
        closest_to_distance_index(&query, &candidates, f64::NAN),
        Err(SpatialError::NonFiniteInput)
    );
}

#[test]
fn average_nearest_neighbour_distances_cover_empty_single_and_duplicates() {
    assert_eq!(average_nearest_neighbour_distance(&[]).to_bits(), 0);
    assert_eq!(
        average_nearest_neighbour_distance(&[point(3.0, 4.0, 5.0)]).to_bits(),
        0
    );
    assert_eq!(
        average_nearest_neighbour_distance(&[
            point(0.0, 0.0, 0.0),
            point(0.0, 0.0, 0.0),
            point(4.0, 0.0, 0.0),
        ]),
        4.0 / 3.0
    );

    let constraints = Constraints {
        inequalities: vec![
            Inequality::new(0.0, 0.0, 0.0, 1.0).unwrap(),
            Inequality::new(4.0, 0.0, 0.0, 2.0).unwrap(),
        ],
        interfaces: vec![
            Interface::new(0.0, 0.0, 0.0, 1.0).unwrap(),
            Interface::new(3.0, 0.0, 0.0, 2.0).unwrap(),
            Interface::new(9.0, 0.0, 0.0, 3.0).unwrap(),
        ],
        planars: vec![
            Planar::from_strike_dip_polarity(7.0, 0.0, 0.0, 0.0, 0.0, Polarity::Upright).unwrap(),
        ],
        tangents: Vec::new(),
    };
    assert_eq!(constraints.compute_inequality_avg_nn_distance(), 4.0);
    assert_eq!(constraints.compute_interface_avg_nn_distance(), 4.0);
    assert_eq!(constraints.compute_planar_avg_nn_distance(), 0.0);
    assert_eq!(constraints.compute_tangent_avg_nn_distance(), 0.0);
    let averages = constraints.compute_avg_nn_distances();
    assert_eq!(averages.inequalities, 4.0);
    assert_eq!(averages.interfaces, 4.0);
    assert_eq!(averages.planars, 0.0);
    assert_eq!(averages.tangents, 0.0);
}

#[test]
fn spatial_metrics_remove_collocated_points_and_compute_xyz_bounds() {
    let points = vec![
        Point::with_c(0.0, 0.0, 0.0, 0.0).unwrap(),
        Point::with_c(0.000_5, 0.0, 0.0, 100.0).unwrap(),
        Point::with_c(2.0, 0.0, 0.0, 0.0).unwrap(),
        Point::with_c(4.0, 0.0, 0.0, 0.0).unwrap(),
    ];
    let metrics = spatial_metrics(&points).unwrap();
    assert_eq!(metrics.resolution, 1.0);
    assert_eq!(metrics.bounds(), [0.0, 4.0, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(bounds(&points).unwrap(), [0.0, 4.0, 0.0, 0.0, 0.0, 0.0]);

    let single = spatial_metrics(&[point(-2.0, 3.0, 7.0)]).unwrap();
    assert_eq!(single.resolution.to_bits(), 0);
    assert_eq!(single.bounds(), [-2.0, -2.0, 3.0, 3.0, 7.0, 7.0]);
    assert_eq!(spatial_metrics(&[]), Err(SpatialError::EmptyPointSet));
    assert_eq!(bounds(&[]), Err(SpatialError::EmptyPointSet));
}

#[test]
fn pair_largest_distance_translation_and_rotation_properties_hold() {
    let points = vec![
        point(0.0, 0.0, 0.0),
        point(3.0, 4.0, 0.0),
        point(-3.0, -4.0, 0.0),
    ];
    assert_eq!(farthest_pair_indices(&points), Ok([1, 2]));
    assert_eq!(largest_distance_between_points(&points), 10.0);
    assert_eq!(
        farthest_pair_indices(&[point(1.0, 1.0, 1.0)]),
        Err(SpatialError::FewerThanTwoPoints)
    );
    assert_eq!(
        farthest_pair_indices(&[]),
        Err(SpatialError::FewerThanTwoPoints)
    );
    assert_eq!(
        farthest_pair_indices(&[point(0.0, 0.0, 0.0), point(0.0, 0.0, 0.0)]),
        Ok([0, 0])
    );
    assert_eq!(largest_distance_between_points(&[]).to_bits(), 0);

    let translated = points
        .iter()
        .map(|value| point(value.x() + 20.0, value.y() - 30.0, value.z() + 7.0))
        .collect::<Vec<_>>();
    let rotated = points
        .iter()
        .map(|value| point(-value.y(), value.x(), value.z()))
        .collect::<Vec<_>>();
    assert_close(
        average_nearest_neighbour_distance(&translated),
        average_nearest_neighbour_distance(&points),
    );
    assert_close(
        largest_distance_between_points(&rotated),
        largest_distance_between_points(&points),
    );
}

#[test]
fn extremal_indices_and_axis_variability_match_frozen_index_order() {
    let points = vec![
        point(0.0, 2.0, 1.0),
        point(10.0, 0.0, 2.0),
        point(5.0, 9.0, -1.0),
        point(3.0, 3.0, 8.0),
        point(7.0, 4.0, 4.0),
        point(4.0, 6.0, 5.0),
        point(6.0, 5.0, 3.0),
    ];
    assert_eq!(extremal_point_indices(&points).unwrap(), [0, 1, 2, 3, 4, 5]);
    assert_eq!(
        extremal_point_indices(&[point(1.0, 1.0, 1.0)]).unwrap(),
        [0]
    );
    assert_eq!(
        extremal_point_indices(&[]),
        Err(SpatialError::EmptyPointSet)
    );

    assert_eq!(
        maximal_axial_variability_order(&[0.0, 10.0, 0.0, 5.0, 0.0, 5.0]).unwrap(),
        [Axis::X, Axis::Z, Axis::Y]
    );
    assert_eq!(
        maximal_axial_variability_order(&[1.0, 1.0, 2.0, 2.0, 3.0, 3.0]).unwrap(),
        [Axis::Z, Axis::Y, Axis::X]
    );
    assert_eq!(
        maximal_axial_variability_order(&[0.0, f64::INFINITY, 0.0, 1.0, 0.0, 1.0]),
        Err(SpatialError::NonFiniteInput)
    );
}
