//! Four-point first-order Lagrangian polynomial basis from frozen Surfe.
//!
//! Sources:
//! - `surfe_lib/basis.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`Lagrangian_Polynomial_Basis`)
//! - `surfe_lib/basis.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`_get_unisolvent_subset`, `_initialize_basis`, `poly`, and `poly_d*`)
//! - `surfe_lib/grbf_exceptions.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`failurecreatinglagrangianpolynomialbasis`)

use crate::ordering::sort_values_with_indices;
use crate::{Error, Interface, Point, POSITION_EPSILON};

#[derive(Clone, Debug)]
struct SelectedPoint {
    point: Point,
    source_index: usize,
}

/// Surfe's four-function, first-order Lagrangian polynomial basis.
///
/// The constructor selects from the first interface group having the strict
/// largest size, then follows the frozen coordinate-range and indexed-sort
/// control flow. The selected horizon and source indices remain available as
/// audit evidence. Axis-aligned planar selections retain Surfe's synthetic
/// `1e-3` adjustment on the first copied point.
#[derive(Clone, Debug)]
pub struct LagrangianPolynomialBasis {
    selected_horizon_index: usize,
    selected_source_indices: [usize; 4],
    unisolvent_points: [Point; 4],
    /// Each row is `[constant, x, y, z]` for one basis function.
    coefficients: [[f64; 4]; 4],
}

impl LagrangianPolynomialBasis {
    /// Select four points and initialize the first-order basis.
    ///
    /// Empty groups, fewer than four points in the largest group, frozen
    /// equal-range selection failures, and zero/non-finite determinants map to
    /// Surfe's Lagrangian-basis construction error. A nonzero determinant is
    /// always attempted, without a condition-number or magnitude threshold.
    pub fn new(interface_point_lists: &[Vec<Interface>]) -> Result<Self, Error> {
        let selected_horizon_index = largest_horizon_index(interface_point_lists)
            .ok_or(Error::LagrangianBasisCreationFailure)?;
        let horizon = &interface_point_lists[selected_horizon_index];
        if horizon.len() < 4 {
            return Err(Error::LagrangianBasisCreationFailure);
        }

        let mut selected = select_subset(horizon)?;
        apply_axis_plane_fallback(horizon, &mut selected)?;

        let selected_source_indices = selected
            .iter()
            .map(|selected| selected.source_index)
            .collect::<Vec<_>>()
            .try_into()
            .map_err(|_| Error::LagrangianBasisCreationFailure)?;
        let unisolvent_points = selected
            .into_iter()
            .map(|selected| selected.point)
            .collect::<Vec<_>>()
            .try_into()
            .map_err(|_| Error::LagrangianBasisCreationFailure)?;
        let coefficients = initialize_coefficients(&unisolvent_points)?;

        Ok(Self {
            selected_horizon_index,
            selected_source_indices,
            unisolvent_points,
            coefficients,
        })
    }

    /// Index of the first interface group with the strict largest point count.
    pub const fn selected_horizon_index(&self) -> usize {
        self.selected_horizon_index
    }

    /// Original indices, within the selected horizon, of the four final points.
    pub const fn selected_source_indices(&self) -> [usize; 4] {
        self.selected_source_indices
    }

    /// Selected point copies used by the basis, including any axis fallback.
    pub const fn unisolvent_points(&self) -> &[Point; 4] {
        &self.unisolvent_points
    }

    /// Basis coefficients in `[constant, x, y, z]` row order.
    pub const fn coefficients(&self) -> &[[f64; 4]; 4] {
        &self.coefficients
    }

    /// Evaluate the four Lagrangian basis functions.
    pub fn values(&self, point: &Point) -> [f64; 4] {
        self.coefficients.map(|coefficient| {
            coefficient[0]
                + coefficient[1] * point.x()
                + coefficient[2] * point.y()
                + coefficient[3] * point.z()
        })
    }

    /// Evaluate derivatives with respect to x in basis-function order.
    pub fn dx(&self, _point: &Point) -> [f64; 4] {
        self.coefficients.map(|coefficient| coefficient[1])
    }

    /// Evaluate derivatives with respect to y in basis-function order.
    pub fn dy(&self, _point: &Point) -> [f64; 4] {
        self.coefficients.map(|coefficient| coefficient[2])
    }

    /// Evaluate derivatives with respect to z in basis-function order.
    pub fn dz(&self, _point: &Point) -> [f64; 4] {
        self.coefficients.map(|coefficient| coefficient[3])
    }
}

fn largest_horizon_index(interface_point_lists: &[Vec<Interface>]) -> Option<usize> {
    let mut selected = (!interface_point_lists.is_empty()).then_some(0)?;
    for index in 1..interface_point_lists.len() {
        if interface_point_lists[index].len() > interface_point_lists[selected].len() {
            selected = index;
        }
    }
    Some(selected)
}

fn select_subset(horizon: &[Interface]) -> Result<Vec<SelectedPoint>, Error> {
    let point_count = horizon.len();
    let mut x_coordinates = horizon
        .iter()
        .map(|value| value.point().x())
        .collect::<Vec<_>>();
    let mut y_coordinates = horizon
        .iter()
        .map(|value| value.point().y())
        .collect::<Vec<_>>();
    let mut z_coordinates = horizon
        .iter()
        .map(|value| value.point().z())
        .collect::<Vec<_>>();
    let mut x_indices = (0..point_count).collect::<Vec<_>>();
    let mut y_indices = x_indices.clone();
    let mut z_indices = x_indices.clone();

    if !sort_values_with_indices(&mut x_coordinates, &mut x_indices)
        || !sort_values_with_indices(&mut y_coordinates, &mut y_indices)
        || !sort_values_with_indices(&mut z_coordinates, &mut z_indices)
    {
        return Err(Error::LagrangianBasisCreationFailure);
    }

    let dx = x_coordinates[point_count - 1] - x_coordinates[0];
    let dy = y_coordinates[point_count - 1] - y_coordinates[0];
    let dz = z_coordinates[point_count - 1] - z_coordinates[0];
    let mut selected = Vec::with_capacity(12);

    // These are deliberately independent conditions. Frozen Surfe therefore
    // selects more than four points and fails when maximum ranges tie.
    if dx >= dy && dx >= dz {
        let secondary = if dy >= dz { &y_indices } else { &z_indices };
        append_axis_selection(horizon, &x_indices, secondary, &mut selected);
    }
    if dy >= dx && dy >= dz {
        let secondary = if dx >= dz { &x_indices } else { &z_indices };
        append_axis_selection(horizon, &y_indices, secondary, &mut selected);
    }
    if dz >= dx && dz >= dy {
        let secondary = if dx >= dy { &x_indices } else { &y_indices };
        append_axis_selection(horizon, &z_indices, secondary, &mut selected);
    }

    if selected.len() == 4 {
        Ok(selected)
    } else {
        Err(Error::LagrangianBasisCreationFailure)
    }
}

fn append_axis_selection(
    horizon: &[Interface],
    primary_indices: &[usize],
    secondary_indices: &[usize],
    selected: &mut Vec<SelectedPoint>,
) {
    let first = primary_indices[0];
    let last = primary_indices[primary_indices.len() - 1];
    push_selected(horizon, first, selected);
    push_selected(horizon, last, selected);

    if let Some(&index) = secondary_indices
        .iter()
        .find(|&&index| index != first && index != last)
    {
        push_selected(horizon, index, selected);
    }
    if let Some(&index) = secondary_indices
        .iter()
        .rev()
        .find(|&&index| index != first && index != last)
    {
        push_selected(horizon, index, selected);
    }
}

fn push_selected(horizon: &[Interface], index: usize, selected: &mut Vec<SelectedPoint>) {
    selected.push(SelectedPoint {
        point: horizon[index].point().clone(),
        source_index: index,
    });
}

fn apply_axis_plane_fallback(
    horizon: &[Interface],
    selected: &mut [SelectedPoint],
) -> Result<(), Error> {
    let first = &selected[0].point;
    let planes = [
        selected[1..]
            .iter()
            .all(|value| value.point.x() == first.x()),
        selected[1..]
            .iter()
            .all(|value| value.point.y() == first.y()),
        selected[1..]
            .iter()
            .all(|value| value.point.z() == first.z()),
    ];
    if !planes.into_iter().any(|on_plane| on_plane) {
        return Ok(());
    }

    let mut found_unique_point = false;
    for (axis, on_plane) in planes.into_iter().enumerate() {
        if !on_plane {
            continue;
        }
        let first_coordinate = coordinate(&selected[0].point, axis);
        if let Some((source_index, replacement)) = horizon
            .iter()
            .enumerate()
            .find(|(_, value)| coordinate(value.point(), axis) != first_coordinate)
        {
            let last = selected
                .last_mut()
                .ok_or(Error::LagrangianBasisCreationFailure)?;
            *last = SelectedPoint {
                point: replacement.point().clone(),
                source_index,
            };
            found_unique_point = true;
        }
    }

    if !found_unique_point {
        let first = &selected[0].point;
        let mut position = first.position();
        for (axis, on_plane) in planes.into_iter().enumerate() {
            if on_plane {
                position[axis] += POSITION_EPSILON;
            }
        }
        selected[0].point = Point::with_c(position[0], position[1], position[2], first.c())
            .map_err(|_| Error::LagrangianBasisCreationFailure)?;
    }

    Ok(())
}

fn coordinate(point: &Point, axis: usize) -> f64 {
    match axis {
        0 => point.x(),
        1 => point.y(),
        2 => point.z(),
        _ => unreachable!("three-dimensional axis index"),
    }
}

#[allow(clippy::many_single_char_names)]
fn initialize_coefficients(points: &[Point; 4]) -> Result<[[f64; 4]; 4], Error> {
    let x1 = points[0].x();
    let y1 = points[0].y();
    let z1 = points[0].z();
    let x2 = points[1].x();
    let y2 = points[1].y();
    let z2 = points[1].z();
    let x3 = points[2].x();
    let y3 = points[2].y();
    let z3 = points[2].z();
    let x4 = points[3].x();
    let y4 = points[3].y();
    let z4 = points[3].z();

    let d = x1 * (y4 * z2 - y3 * z2 + y2 * z3 - y4 * z3 - y2 * z4 + y3 * z4)
        + x2 * (y3 * z1 - y4 * z1 - y1 * z3 + y4 * z3 + y1 * z4 - y3 * z4)
        + x3 * (y4 * z1 + y1 * z2 - y4 * z2 - y1 * z4 - y2 * z1 + y2 * z4)
        + x4 * (y2 * z1 - y3 * z1 - y1 * z2 + y3 * z2 + y1 * z3 - y2 * z3);
    if d == 0.0 || !d.is_finite() {
        return Err(Error::LagrangianBasisCreationFailure);
    }

    let coefficients = [
        [
            (x4 * y3 * z2 - x3 * y4 * z2 - x4 * y2 * z3 + x2 * y4 * z3 + x3 * y2 * z4
                - x2 * y3 * z4)
                / d,
            (-(y3 * z2) + y4 * z2 + y2 * z3 - y4 * z3 - y2 * z4 + y3 * z4) / d,
            (x3 * z2 - x4 * z2 - x2 * z3 + x4 * z3 + x2 * z4 - x3 * z4) / d,
            (-(x3 * y2) + x4 * y2 + x2 * y3 - x4 * y3 - x2 * y4 + x3 * y4) / d,
        ],
        [
            (-(x4 * y3 * z1) + x3 * y4 * z1 + x4 * y1 * z3 - x1 * y4 * z3 - x3 * y1 * z4
                + x1 * y3 * z4)
                / d,
            (y3 * z1 - y4 * z1 - y1 * z3 + y4 * z3 + y1 * z4 - y3 * z4) / d,
            (-(x3 * z1) + x4 * z1 + x1 * z3 - x4 * z3 - x1 * z4 + x3 * z4) / d,
            (x3 * y1 - x4 * y1 - x1 * y3 + x4 * y3 + x1 * y4 - x3 * y4) / d,
        ],
        [
            (x4 * y2 * z1 - x2 * y4 * z1 - x4 * y1 * z2 + x1 * y4 * z2 + x2 * y1 * z4
                - x1 * y2 * z4)
                / d,
            (-(y2 * z1) + y4 * z1 + y1 * z2 - y4 * z2 - y1 * z4 + y2 * z4) / d,
            (x2 * z1 - x4 * z1 - x1 * z2 + x4 * z2 + x1 * z4 - x2 * z4) / d,
            (-(x2 * y1) + x4 * y1 + x1 * y2 - x4 * y2 - x1 * y4 + x2 * y4) / d,
        ],
        [
            (-(x3 * y2 * z1) + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3
                + x1 * y2 * z3)
                / d,
            (y2 * z1 - y3 * z1 - y1 * z2 + y3 * z2 + y1 * z3 - y2 * z3) / d,
            (-(x2 * z1) + x3 * z1 + x1 * z2 - x3 * z2 - x1 * z3 + x2 * z3) / d,
            (x2 * y1 - x3 * y1 - x1 * y2 + x3 * y2 + x1 * y3 - x2 * y3) / d,
        ],
    ];

    if coefficients.iter().flatten().all(|value| value.is_finite()) {
        Ok(coefficients)
    } else {
        Err(Error::LagrangianBasisCreationFailure)
    }
}
