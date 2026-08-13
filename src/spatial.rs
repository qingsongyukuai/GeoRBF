//! Surfe-compatible spatial helper algorithms.
//!
//! Source: `surfe_lib/modelling_input.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`convert_constraints_to_points`, `distance_btw_pts`, nearest/farthest
//! helpers, `avg_nn_distance`, `spatial_metrics`, bounds/extremal helpers, and
//! `get_maximal_axial_variability_order`).

use std::fmt;

use crate::{collocated, compare_points, sort_values_with_indices, Axis, Constraints, Point};

/// Safe failures for spatial inputs that frozen Surfe represents with an
/// invalid index, sentinel bounds, or an out-of-bounds access.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum SpatialError {
    /// A point slice required by the operation is empty.
    EmptyPointSet,
    /// A pair operation received fewer than two points.
    FewerThanTwoPoints,
    /// No candidate has a non-zero distance from the query point.
    NoNonzeroNeighbour,
    /// A scalar or explicit bounds value is NaN or infinite.
    NonFiniteInput,
    /// Frozen Surfe's indexed sort could not complete.
    IndexedSortFailure,
}

impl fmt::Display for SpatialError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::EmptyPointSet => "spatial operation requires at least one point",
            Self::FewerThanTwoPoints => "spatial pair operation requires at least two points",
            Self::NoNonzeroNeighbour => "no non-zero-distance neighbour exists",
            Self::NonFiniteInput => "spatial input must be finite",
            Self::IndexedSortFailure => "indexed spatial sort failed",
        })
    }
}

impl std::error::Error for SpatialError {}

/// Bounds and half-average-nearest-neighbour resolution returned by Surfe's
/// `spatial_metrics`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialParameters {
    pub resolution: f64,
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub zmin: f64,
    pub zmax: f64,
}

impl SpatialParameters {
    /// Return `[xmin, xmax, ymin, ymax, zmin, zmax]`.
    pub const fn bounds(self) -> [f64; 6] {
        [
            self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax,
        ]
    }
}

/// Per-category values computed by `Constraints::compute_avg_nn_distances`.
///
/// Frozen Surfe caches these four values in mutable constraint state. GeoRBF
/// returns an owned snapshot, keeping the spatial calculation pure and free of
/// hidden shared mutation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConstraintAverageNearestNeighbourDistances {
    pub inequalities: f64,
    pub interfaces: f64,
    pub planars: f64,
    pub tangents: f64,
}

/// Convert the four constraint categories to points in frozen category order:
/// inequality, interface, planar, tangent.
pub fn constraints_to_points(constraints: &Constraints) -> Vec<Point> {
    constraints
        .inequalities
        .iter()
        .map(|value| value.point().clone())
        .chain(
            constraints
                .interfaces
                .iter()
                .map(|value| value.point().clone()),
        )
        .chain(
            constraints
                .planars
                .iter()
                .map(|value| value.point().clone()),
        )
        .chain(
            constraints
                .tangents
                .iter()
                .map(|value| value.point().clone()),
        )
        .collect()
}

/// Four-dimensional Euclidean distance, including Surfe's `c` coordinate.
pub fn distance_between_points(first: &Point, second: &Point) -> f64 {
    let dx = first.x() - second.x();
    let dy = first.y() - second.y();
    let dz = first.z() - second.z();
    let dc = first.c() - second.c();
    (dx * dx + dy * dy + dz * dz + dc * dc).sqrt()
}

/// Index of the first nearest candidate with non-zero distance.
pub fn nearest_neighbour_index(query: &Point, points: &[Point]) -> Result<usize, SpatialError> {
    if points.is_empty() {
        return Err(SpatialError::EmptyPointSet);
    }

    let mut minimum_distance = f64::MAX;
    let mut nearest = None;
    for (index, point) in points.iter().enumerate() {
        let distance = distance_between_points(query, point);
        if distance != 0.0 && distance < minimum_distance {
            minimum_distance = distance;
            nearest = Some(index);
        }
    }
    nearest.ok_or(SpatialError::NoNonzeroNeighbour)
}

/// Up to `count` nearest non-zero-distance candidate indices.
///
/// Distances use frozen `sort_vector_w_index`, including its tie ordering.
/// Frozen C++ can index past its filtered arrays when `count` includes exact
/// query duplicates; the safe Rust result is capped to the candidates that
/// actually remain.
pub fn nearest_neighbour_indices(count: usize, query: &Point, points: &[Point]) -> Vec<usize> {
    let mut distances = Vec::with_capacity(points.len());
    let mut indices = Vec::with_capacity(points.len());
    for (index, point) in points.iter().enumerate() {
        let distance = distance_between_points(query, point);
        if distance != 0.0 {
            distances.push(distance);
            indices.push(index);
        }
    }
    let sorted = sort_values_with_indices(&mut distances, &mut indices);
    debug_assert!(
        sorted,
        "parallel distance and index vectors have equal lengths"
    );
    indices.truncate(count.min(indices.len()));
    indices
}

/// Index of the first farthest point from a query.
pub fn farthest_neighbour_index(query: &Point, points: &[Point]) -> Result<usize, SpatialError> {
    let first = points.first().ok_or(SpatialError::EmptyPointSet)?;
    let mut index = 0;
    let mut largest_distance = distance_between_points(query, first);
    for (candidate_index, point) in points.iter().enumerate().skip(1) {
        let distance = distance_between_points(query, point);
        if distance > largest_distance {
            largest_distance = distance;
            index = candidate_index;
        }
    }
    Ok(index)
}

/// Index in `points` having the largest pairwise distance to any member of
/// `other_points`.
pub fn farthest_from_other_set_index(
    points: &[Point],
    other_points: &[Point],
) -> Result<usize, SpatialError> {
    if points.is_empty() || other_points.is_empty() {
        return Err(SpatialError::EmptyPointSet);
    }

    let mut index = 0;
    let mut largest_distance = 0.0;
    for (candidate_index, point) in points.iter().enumerate() {
        for other in other_points {
            let distance = distance_between_points(point, other);
            if distance > largest_distance {
                largest_distance = distance;
                index = candidate_index;
            }
        }
    }
    Ok(index)
}

/// Mean distance from every point to its nearest *other-index* point.
///
/// Empty and singleton inputs both return zero. Exact duplicate points have a
/// nearest distance of zero, matching the source's `k != j` rule.
pub fn average_nearest_neighbour_distance(points: &[Point]) -> f64 {
    let mut average = 0.0;
    let count = points.len();
    for (index, point) in points.iter().enumerate() {
        let mut minimum_distance = f64::MAX;
        for (other_index, other) in points.iter().enumerate() {
            if other_index != index {
                let distance = distance_between_points(point, other);
                if distance < minimum_distance {
                    minimum_distance = distance;
                }
            }
        }
        if count == 1 {
            minimum_distance = 0.0;
        }
        average += minimum_distance;
    }
    if count != 0 {
        average /= count as f64;
    }
    average
}

/// Remove collocated positions, compute XYZ bounds, and return half the
/// average nearest-neighbour distance as resolution.
pub fn spatial_metrics(points: &[Point]) -> Result<SpatialParameters, SpatialError> {
    if points.is_empty() {
        return Err(SpatialError::EmptyPointSet);
    }

    let mut distinct = points.to_vec();
    distinct.sort_by(compare_points);
    distinct.dedup_by(|left, right| collocated(left, right));
    let point_bounds = bounds(&distinct)?;
    Ok(SpatialParameters {
        resolution: average_nearest_neighbour_distance(&distinct) / 2.0,
        xmin: point_bounds[0],
        xmax: point_bounds[1],
        ymin: point_bounds[2],
        ymax: point_bounds[3],
        zmin: point_bounds[4],
        zmax: point_bounds[5],
    })
}

/// Ordered indices of the first pair with the largest distance.
pub fn farthest_pair_indices(points: &[Point]) -> Result<[usize; 2], SpatialError> {
    if points.len() < 2 {
        return Err(SpatialError::FewerThanTwoPoints);
    }

    let mut indices = [0, 0];
    let mut largest_distance = -f64::MAX;
    for (first_index, first) in points.iter().enumerate() {
        for (second_index, second) in points.iter().enumerate() {
            let distance = distance_between_points(first, second);
            if distance > largest_distance {
                largest_distance = distance;
                indices = [first_index, second_index];
            }
        }
    }
    Ok(indices)
}

/// Index whose query distance is first-closest to `target_distance`.
pub fn closest_to_distance_index(
    query: &Point,
    points: &[Point],
    target_distance: f64,
) -> Result<usize, SpatialError> {
    if points.is_empty() {
        return Err(SpatialError::EmptyPointSet);
    }
    if !target_distance.is_finite() {
        return Err(SpatialError::NonFiniteInput);
    }

    let mut smallest_residual = f64::MAX;
    let mut closest = None;
    for (index, point) in points.iter().enumerate() {
        let residual = (distance_between_points(query, point) - target_distance).abs();
        if residual < smallest_residual {
            smallest_residual = residual;
            closest = Some(index);
        }
    }
    closest.ok_or(SpatialError::NonFiniteInput)
}

/// XYZ bounds in frozen `[xmin, xmax, ymin, ymax, zmin, zmax]` order.
pub fn bounds(points: &[Point]) -> Result<[f64; 6], SpatialError> {
    let first = points.first().ok_or(SpatialError::EmptyPointSet)?;
    let mut result = [
        first.x(),
        first.x(),
        first.y(),
        first.y(),
        first.z(),
        first.z(),
    ];
    for point in &points[1..] {
        if point.x() < result[0] {
            result[0] = point.x();
        }
        if point.x() > result[1] {
            result[1] = point.x();
        }
        if point.y() < result[2] {
            result[2] = point.y();
        }
        if point.y() > result[3] {
            result[3] = point.y();
        }
        if point.z() < result[4] {
            result[4] = point.z();
        }
        if point.z() > result[5] {
            result[5] = point.z();
        }
    }
    Ok(result)
}

/// Indices that sample coordinate extrema in descending axial-range order.
pub fn extremal_point_indices(points: &[Point]) -> Result<Vec<usize>, SpatialError> {
    if points.is_empty() {
        return Err(SpatialError::EmptyPointSet);
    }

    let count = points.len();
    let mut x_coordinates = points.iter().map(Point::x).collect::<Vec<_>>();
    let mut y_coordinates = points.iter().map(Point::y).collect::<Vec<_>>();
    let mut z_coordinates = points.iter().map(Point::z).collect::<Vec<_>>();
    let mut x_indices = (0..count).collect::<Vec<_>>();
    let mut y_indices = x_indices.clone();
    let mut z_indices = x_indices.clone();
    if !sort_values_with_indices(&mut x_coordinates, &mut x_indices)
        || !sort_values_with_indices(&mut y_coordinates, &mut y_indices)
        || !sort_values_with_indices(&mut z_coordinates, &mut z_indices)
    {
        return Err(SpatialError::IndexedSortFailure);
    }

    let mut ranges = vec![
        x_coordinates[count - 1] - x_coordinates[0],
        y_coordinates[count - 1] - y_coordinates[0],
        z_coordinates[count - 1] - z_coordinates[0],
    ];
    let mut axes = vec![0_usize, 1, 2];
    if !sort_values_with_indices(&mut ranges, &mut axes) {
        return Err(SpatialError::IndexedSortFailure);
    }

    let mut result = Vec::with_capacity(count.min(6));
    for axis in axes.into_iter().rev() {
        let sorted_indices = match axis {
            0 => &x_indices,
            1 => &y_indices,
            2 => &z_indices,
            _ => unreachable!("axis indices are initialized to 0, 1, 2"),
        };
        append_unused_extreme(sorted_indices, false, &mut result);
        append_unused_extreme(sorted_indices, true, &mut result);
    }
    Ok(result)
}

fn append_unused_extreme(sorted: &[usize], from_high_end: bool, result: &mut Vec<usize>) {
    let candidate = if from_high_end {
        sorted
            .iter()
            .rev()
            .copied()
            .find(|index| !result.contains(index))
    } else {
        sorted.iter().copied().find(|index| !result.contains(index))
    };
    if let Some(index) = candidate {
        result.push(index);
    }
}

/// Axes ordered from largest to smallest absolute bounds range.
pub fn maximal_axial_variability_order(bounds: &[f64; 6]) -> Result<[Axis; 3], SpatialError> {
    if !bounds.iter().all(|value| value.is_finite()) {
        return Err(SpatialError::NonFiniteInput);
    }

    let mut ranges = vec![
        (bounds[0] - bounds[1]).abs(),
        (bounds[2] - bounds[3]).abs(),
        (bounds[4] - bounds[5]).abs(),
    ];
    let mut axes = vec![0_usize, 1, 2];
    if !sort_values_with_indices(&mut ranges, &mut axes) {
        return Err(SpatialError::IndexedSortFailure);
    }
    Ok([
        axis_from_index(axes[2]),
        axis_from_index(axes[1]),
        axis_from_index(axes[0]),
    ])
}

fn axis_from_index(index: usize) -> Axis {
    match index {
        0 => Axis::X,
        1 => Axis::Y,
        2 => Axis::Z,
        _ => unreachable!("axis indices are initialized to 0, 1, 2"),
    }
}

/// Largest distance over all ordered point pairs; zero for fewer than two
/// distinct positions, including empty and singleton inputs.
pub fn largest_distance_between_points(points: &[Point]) -> f64 {
    let mut largest_distance = 0.0;
    for first in points {
        for second in points {
            let distance = distance_between_points(first, second);
            if distance > largest_distance {
                largest_distance = distance;
            }
        }
    }
    largest_distance
}

impl Constraints {
    pub fn compute_inequality_avg_nn_distance(&self) -> f64 {
        let points = self
            .inequalities
            .iter()
            .map(|value| value.point().clone())
            .collect::<Vec<_>>();
        average_nearest_neighbour_distance(&points)
    }

    pub fn compute_interface_avg_nn_distance(&self) -> f64 {
        let points = self
            .interfaces
            .iter()
            .map(|value| value.point().clone())
            .collect::<Vec<_>>();
        average_nearest_neighbour_distance(&points)
    }

    pub fn compute_planar_avg_nn_distance(&self) -> f64 {
        let points = self
            .planars
            .iter()
            .map(|value| value.point().clone())
            .collect::<Vec<_>>();
        average_nearest_neighbour_distance(&points)
    }

    pub fn compute_tangent_avg_nn_distance(&self) -> f64 {
        let points = self
            .tangents
            .iter()
            .map(|value| value.point().clone())
            .collect::<Vec<_>>();
        average_nearest_neighbour_distance(&points)
    }

    pub fn compute_avg_nn_distances(&self) -> ConstraintAverageNearestNeighbourDistances {
        ConstraintAverageNearestNeighbourDistances {
            inequalities: self.compute_inequality_avg_nn_distance(),
            interfaces: self.compute_interface_avg_nn_distance(),
            planars: self.compute_planar_avg_nn_distance(),
            tangents: self.compute_tangent_avg_nn_distance(),
        }
    }
}
