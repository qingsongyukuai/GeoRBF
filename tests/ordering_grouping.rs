use std::cmp::Ordering;

use georbf::{
    collocated, compare_points, sort_values_with_indices, Constraints, Inequality, Interface,
    Planar, Point, Tangent,
};

fn interface(x: f64, y: f64, z: f64, level: f64) -> Interface {
    Interface::new(x, y, z, level).unwrap()
}

fn grouping_fixture() -> Constraints {
    Constraints {
        interfaces: vec![
            interface(0.08, 0.0, 0.0, -0.0),
            interface(0.02, 0.0, 0.0, 2.0),
            interface(0.07, 0.0, 0.0, 0.0),
            interface(0.05, 0.0, 0.0, f64::from_bits(1.0_f64.to_bits() + 1)),
            interface(0.03, 0.0, 0.0, 2.0),
            interface(0.01, 0.0, 0.0, 3.0),
            interface(0.06, 0.0, 0.0, 1.0),
            interface(0.04, 0.0, 0.0, 1.0),
        ],
        ..Constraints::default()
    }
}

fn grouped_positions(constraints: &Constraints) -> (Vec<u64>, Vec<f64>, Vec<Vec<f64>>) {
    let grouping = constraints.interface_grouping().unwrap();
    let levels = grouping
        .levels_descending()
        .iter()
        .map(|level| level.to_bits())
        .collect();
    let references = grouping
        .reference_indices()
        .iter()
        .map(|&index| constraints.interfaces[index].point().x())
        .collect();
    let groups = grouping
        .multi_point_groups()
        .iter()
        .map(|indices| {
            indices
                .iter()
                .map(|&index| constraints.interfaces[index].point().x())
                .collect()
        })
        .collect();
    (levels, references, groups)
}

#[test]
fn point_order_and_collocation_match_frozen_surfe_boundaries() {
    let origin = Point::new(0.0, -0.0, 0.0).unwrap();
    let below = f64::from_bits(0.001_f64.to_bits() - 1);
    let above = f64::from_bits(0.001_f64.to_bits() + 1);
    assert!(collocated(
        &origin,
        &Point::new(below, below, below).unwrap()
    ));
    assert!(!collocated(
        &origin,
        &Point::new(0.001, 0.001, 0.001).unwrap()
    ));
    assert!(!collocated(
        &origin,
        &Point::new(above, above, above).unwrap()
    ));
    assert!(collocated(
        &Point::with_c(1.0, 2.0, 3.0, -99.0).unwrap(),
        &Point::with_c(1.0, 2.0, 3.0, 99.0).unwrap()
    ));

    assert_eq!(
        compare_points(
            &Point::new(-0.0, 1.0, 0.0).unwrap(),
            &Point::new(0.0, 1.0, -0.0).unwrap()
        ),
        Ordering::Equal
    );
    assert_eq!(
        compare_points(
            &Point::new(0.0, -1.0, 9.0).unwrap(),
            &Point::new(0.0, 1.0, -9.0).unwrap()
        ),
        Ordering::Less
    );

    let mut points = [
        Point::with_c(0.0, -1.0, 0.0, 11.0).unwrap(),
        Point::new(1.0, 0.0, 0.0).unwrap(),
        Point::new(0.0, -1.0, -1.0).unwrap(),
        Point::with_c(-0.0, -1.0, -0.0, 22.0).unwrap(),
        Point::new(-1.0, 5.0, 5.0).unwrap(),
    ];
    points.sort_by(compare_points);
    assert_eq!(
        points.iter().map(Point::position).collect::<Vec<_>>(),
        vec![
            [-1.0, 5.0, 5.0],
            [0.0, -1.0, -1.0],
            [0.0, -1.0, 0.0],
            [-0.0, -1.0, -0.0],
            [1.0, 0.0, 0.0],
        ]
    );
    assert_eq!((points[2].c(), points[3].c()), (11.0, 22.0));
}

#[test]
fn indexed_sort_matches_frozen_numerical_recipes_order() {
    let mut small = vec![3.0, -0.0, 0.0, 2.0, -1.0, 2.0];
    let mut small_indices = (0..small.len()).collect::<Vec<_>>();
    assert!(sort_values_with_indices(&mut small, &mut small_indices));
    assert_eq!(
        small
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        [-1.0_f64, -0.0, 0.0, 2.0, 2.0, 3.0]
            .into_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>()
    );
    assert_eq!(small_indices, [4, 1, 2, 3, 5, 0]);

    let mut partition = vec![5.0, 1.0, 3.0, 1.0, 2.0, 5.0, 0.0, -0.0, 4.0, 3.0, -1.0, 2.0];
    let mut partition_indices = (0..partition.len()).collect::<Vec<_>>();
    assert!(sort_values_with_indices(
        &mut partition,
        &mut partition_indices
    ));
    assert_eq!(partition_indices, [10, 6, 7, 3, 1, 11, 4, 2, 9, 8, 5, 0]);
    assert_eq!(partition[1].to_bits(), 0.0_f64.to_bits());
    assert_eq!(partition[2].to_bits(), (-0.0_f64).to_bits());

    let mut mismatch_values = [2.0, 1.0];
    let mut mismatch_indices = [9];
    assert!(!sort_values_with_indices(
        &mut mismatch_values,
        &mut mismatch_indices
    ));
    assert_eq!(mismatch_values, [2.0, 1.0]);
    assert_eq!(mismatch_indices, [9]);
}

#[test]
fn removal_is_sorted_category_local_and_keeps_the_first_survivor() {
    let below = f64::from_bits(0.001_f64.to_bits() - 1);
    let above = f64::from_bits(0.001_f64.to_bits() + 1);
    let mut constraints = Constraints {
        interfaces: vec![
            interface(above, 0.0, 0.0, 444.0),
            interface(0.001, 0.0, 0.0, 333.0),
            interface(below, 0.0, 0.0, 222.0),
            interface(0.0, 0.0, 0.0, 111.0),
        ],
        inequalities: vec![
            Inequality::new(9.0, 9.0, 9.0, 10.0).unwrap(),
            Inequality::new(9.0, 9.0, 9.0, 20.0).unwrap(),
        ],
        planars: vec![
            Planar::from_normal(9.0, 9.0, 9.0, 1.0, 0.0, 0.0).unwrap(),
            Planar::from_normal(9.0, 9.0, 9.0, 0.0, 1.0, 0.0).unwrap(),
        ],
        tangents: vec![
            Tangent::new(9.0, 9.0, 9.0, 1.0, 0.0, 0.0).unwrap(),
            Tangent::new(9.0, 9.0, 9.0, 0.0, 1.0, 0.0).unwrap(),
        ],
    };

    let removed = constraints.remove_collocated();
    assert_eq!(removed.interfaces, 2);
    assert_eq!(removed.inequalities, 1);
    assert_eq!(removed.planars, 1);
    assert_eq!(removed.tangents, 1);
    assert_eq!(
        constraints
            .interfaces
            .iter()
            .map(|value| (value.point().x().to_bits(), value.level()))
            .collect::<Vec<_>>(),
        vec![(0.0_f64.to_bits(), 111.0), (0.001_f64.to_bits(), 333.0)]
    );
    assert_eq!(constraints.inequalities[0].level(), 10.0);
    assert_eq!(constraints.planars[0].normal(), [1.0, 0.0, 0.0]);
    assert_eq!(constraints.tangents[0].vector(), [1.0, 0.0, 0.0]);

    assert_eq!(constraints.inequalities.len(), 1);
    assert_eq!(constraints.interfaces.len(), 2);
    assert_eq!(constraints.planars.len(), 1);
    assert_eq!(constraints.tangents.len(), 1);

    let mut cross_category = Constraints {
        inequalities: vec![Inequality::new(5.0, 5.0, 5.0, 1.0).unwrap()],
        interfaces: vec![interface(5.0, 5.0, 5.0, 1.0)],
        planars: vec![Planar::from_normal(5.0, 5.0, 5.0, 0.0, 0.0, 1.0).unwrap()],
        tangents: vec![Tangent::new(5.0, 5.0, 5.0, 1.0, 0.0, 0.0).unwrap()],
    };
    assert_eq!(cross_category.remove_collocated().total(), 0);
    assert_eq!(
        (
            cross_category.inequalities.len(),
            cross_category.interfaces.len(),
            cross_category.planars.len(),
            cross_category.tangents.len(),
        ),
        (1, 1, 1, 1)
    );
}

#[test]
fn exact_levels_references_groups_and_counts_match_frozen_surfe() {
    let mut constraints = grouping_fixture();
    assert_eq!(constraints.remove_collocated().total(), 0);

    let grouping = constraints.interface_grouping().unwrap();
    assert_eq!(
        grouping
            .levels_descending()
            .iter()
            .map(|level| level.to_bits())
            .collect::<Vec<_>>(),
        [
            3.0_f64.to_bits(),
            2.0_f64.to_bits(),
            f64::from_bits(1.0_f64.to_bits() + 1).to_bits(),
            1.0_f64.to_bits(),
            0.0_f64.to_bits(),
        ]
    );
    assert_eq!(grouping.reference_indices(), [0, 1, 4, 3, 6]);
    assert_eq!(
        grouping.multi_point_groups(),
        [vec![1, 2], vec![3, 5], vec![6, 7]]
    );
    assert_eq!(grouping.increment_pair_count(), 3);
    assert_eq!(grouping.sequenced_reference_pair_count(), 4);

    constraints.inequalities = vec![
        Inequality::new(0.0, 0.0, 0.0, -0.0).unwrap(),
        Inequality::new(1.0, 0.0, 0.0, 2.0).unwrap(),
        Inequality::new(2.0, 0.0, 0.0, 1.0).unwrap(),
        Inequality::new(3.0, 0.0, 0.0, 2.0).unwrap(),
        Inequality::new(4.0, 0.0, 0.0, 0.0).unwrap(),
    ];
    assert_eq!(
        constraints
            .distinct_inequality_levels()
            .iter()
            .map(|level| level.to_bits())
            .collect::<Vec<_>>(),
        [2.0_f64.to_bits(), 1.0_f64.to_bits(), (-0.0_f64).to_bits()]
    );
}

#[test]
fn grouping_is_reproducible_for_permuted_input_and_empty_input_is_explicit() {
    let mut forward = grouping_fixture();
    let mut reversed = grouping_fixture();
    reversed.interfaces.reverse();
    forward.remove_collocated();
    reversed.remove_collocated();

    assert_eq!(grouped_positions(&forward), grouped_positions(&reversed));

    let empty = Constraints::default();
    assert!(empty.interface_grouping().is_none());
    assert!(empty.distinct_inequality_levels().is_empty());
}
