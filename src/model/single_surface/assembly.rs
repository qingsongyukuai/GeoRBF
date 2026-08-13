//! Frozen Single Surface matrix partitions and right-hand sides.
//!
//! Source: `surfe_lib/single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`get_equality_values`, both `get_inequality_values` overloads, and
//! `get_inequality_matrix`).

use crate::{
    assembly::{normal_component, rows_for_range},
    constraints_to_points, largest_distance_between_points, AssemblyConstraints, AssemblyError,
    BoundedConstraintSystem, ConstraintLayout, ConstraintSystem, Constraints, DenseMatrix,
    DenseVector, Error, LayoutDof, Parameters,
};

pub(crate) fn build(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    interpolation: &DenseMatrix,
) -> Result<AssemblyConstraints, AssemblyError> {
    if layout.internal_parameters().restricted_range {
        return bounded(layout, constraints, parameters, interpolation);
    }
    if layout.internal_parameters().modified_basis {
        return quadratic(layout, constraints, interpolation);
    }
    Ok(AssemblyConstraints::Linear {
        right_hand_side: DenseVector::from_values(
            layout
                .dofs()
                .iter()
                .map(|dof| linear_value(dof, constraints))
                .collect::<Result<Vec<_>, _>>()?,
        ),
    })
}

fn linear_value(dof: &LayoutDof, constraints: &Constraints) -> Result<f64, AssemblyError> {
    match dof {
        LayoutDof::InterfaceValue { index } => constraints
            .interfaces
            .get(*index)
            .map(|value| value.level())
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        LayoutDof::PlanarDerivative { index, axis } => constraints
            .planars
            .get(*index)
            .map(|value| normal_component(value.normal(), *axis))
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        LayoutDof::Tangent { index } => constraints
            .tangents
            .get(*index)
            .map(|value| value.inner_product_constraint())
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        LayoutDof::PolynomialTerm { .. } => Ok(0.0),
        _ => Err(AssemblyError::Surfe(Error::EqualityVectorFailure)),
    }
}

fn quadratic(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    interpolation: &DenseMatrix,
) -> Result<AssemblyConstraints, AssemblyError> {
    let equality_range = layout.partitions().equality();
    let inequality_range = layout.partitions().inequality();
    if equality_range.is_empty() {
        return Err(AssemblyError::Surfe(Error::EqualityVectorFailure));
    }
    if inequality_range.is_empty() {
        return Err(AssemblyError::Surfe(Error::InequalityVectorFailure));
    }
    let equality = ConstraintSystem::new(
        rows_for_range(interpolation, equality_range),
        DenseVector::from_values(
            layout.dofs()[equality_range.start()..equality_range.end()]
                .iter()
                .map(|dof| linear_value(dof, constraints))
                .collect::<Result<Vec<_>, _>>()?,
        ),
    );
    let mut inequality_matrix = rows_for_range(interpolation, inequality_range);
    for (row, dof) in layout.dofs()[inequality_range.start()..inequality_range.end()]
        .iter()
        .enumerate()
    {
        let LayoutDof::InequalityValue { index } = dof else {
            return Err(AssemblyError::Surfe(Error::InequalityVectorFailure));
        };
        let level = constraints
            .inequalities
            .get(*index)
            .ok_or(AssemblyError::Surfe(Error::InequalityVectorFailure))?
            .level();
        if level <= 0.0 {
            for column in 0..inequality_matrix.cols() {
                let value = inequality_matrix
                    .get(row, column)
                    .expect("partition matrix indices are valid");
                inequality_matrix.set(row, column, -value);
            }
        }
    }
    let inequality = ConstraintSystem::new(
        inequality_matrix,
        DenseVector::zeros(inequality_range.len()),
    );
    Ok(AssemblyConstraints::Quadratic {
        equality,
        inequality,
    })
}

fn bounded(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    interpolation: &DenseMatrix,
) -> Result<AssemblyConstraints, AssemblyError> {
    let mut bounded_constraints = constraints.clone();
    for interface in &mut bounded_constraints.interfaces {
        interface.set_level_bounds(parameters.interface_uncertainty)?;
    }
    for planar in &mut bounded_constraints.planars {
        planar.set_normal_bounds(
            parameters.angular_uncertainty,
            parameters.angular_uncertainty / 2.0,
        )?;
    }
    for tangent in &mut bounded_constraints.tangents {
        tangent.set_angle_bounds(parameters.angular_uncertainty)?;
    }
    let points = constraints_to_points(&bounded_constraints);
    let distance = largest_distance_between_points(&points);
    let mut lower = Vec::with_capacity(layout.matrix_size());
    let mut range = Vec::with_capacity(layout.matrix_size());
    for dof in layout.dofs() {
        let bounds = match dof {
            LayoutDof::InequalityValue { index } => {
                let level = bounded_constraints
                    .inequalities
                    .get(*index)
                    .ok_or(AssemblyError::Surfe(Error::InequalityVectorFailure))?
                    .level();
                if level > 0.0 {
                    (0.0, distance)
                } else {
                    (-distance, distance)
                }
            }
            LayoutDof::InterfaceValue { index } => {
                let interface = bounded_constraints
                    .interfaces
                    .get(*index)
                    .ok_or(AssemblyError::Surfe(Error::InequalityVectorFailure))?;
                (
                    interface.level_lower_bound(),
                    interface.level_upper_bound() - interface.level_lower_bound(),
                )
            }
            LayoutDof::PlanarDerivative { index, axis } => {
                let bounds = bounded_constraints
                    .planars
                    .get(*index)
                    .and_then(|value| value.normal_bounds())
                    .ok_or(AssemblyError::Surfe(Error::InequalityVectorFailure))?;
                let axis = match axis {
                    crate::Axis::X => 0,
                    crate::Axis::Y => 1,
                    crate::Axis::Z => 2,
                };
                (bounds[axis][0], bounds[axis][1] - bounds[axis][0])
            }
            LayoutDof::Tangent { index } => {
                let bounds = bounded_constraints
                    .tangents
                    .get(*index)
                    .and_then(|value| value.angle_bounds())
                    .ok_or(AssemblyError::Surfe(Error::InequalityVectorFailure))?;
                (bounds[0], bounds[1] - bounds[0])
            }
            _ => return Err(AssemblyError::Surfe(Error::InequalityVectorFailure)),
        };
        lower.push(bounds.0);
        range.push(bounds.1);
    }
    Ok(AssemblyConstraints::Bounded {
        system: BoundedConstraintSystem::new(
            interpolation.clone(),
            DenseVector::from_values(lower),
            DenseVector::from_values(range),
        ),
    })
}
