//! Frozen Stratigraphic Horizons QP partitions, values, and bounds.
//!
//! Source:
//! `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`get_equality_values`, both `get_inequality_values` overloads, and
//! `get_inequality_matrix`).

use crate::{
    assembly::{normal_component, rows_for_range},
    constraints_to_points, largest_distance_between_points, AssemblyConstraints, AssemblyError,
    BoundedConstraintSystem, ConstraintLayout, ConstraintSystem, Constraints, DenseMatrix,
    DenseVector, DifferenceKind, Error, LayoutDof, LayoutPointRef, Parameters,
};

pub(crate) fn build(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    interpolation: &DenseMatrix,
) -> Result<AssemblyConstraints, AssemblyError> {
    if layout.internal_parameters().restricted_range {
        bounded(layout, constraints, parameters, interpolation)
    } else {
        quadratic(layout, constraints, parameters, interpolation)
    }
}

fn quadratic(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
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
    let equality_values = layout.dofs()[equality_range.start()..equality_range.end()]
        .iter()
        .map(|dof| equality_value(dof, constraints))
        .collect::<Result<Vec<_>, _>>()?;
    let inequality_values = layout.dofs()[inequality_range.start()..inequality_range.end()]
        .iter()
        .map(|dof| inequality_value(dof, parameters))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(AssemblyConstraints::Quadratic {
        equality: ConstraintSystem::new(
            rows_for_range(interpolation, equality_range),
            DenseVector::from_values(equality_values),
        ),
        inequality: ConstraintSystem::new(
            rows_for_range(interpolation, inequality_range),
            DenseVector::from_values(inequality_values),
        ),
    })
}

fn equality_value(dof: &LayoutDof, constraints: &Constraints) -> Result<f64, AssemblyError> {
    match dof {
        LayoutDof::Difference {
            kind: DifferenceKind::SameLevelInterface,
            positive,
            negative,
        } => Ok(level(*positive, constraints)? - level(*negative, constraints)?),
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
        _ => Err(AssemblyError::Surfe(Error::EqualityVectorFailure)),
    }
}

fn inequality_value(dof: &LayoutDof, parameters: &Parameters) -> Result<f64, AssemblyError> {
    match dof {
        LayoutDof::Difference {
            kind: DifferenceKind::SequencedInterfaces,
            ..
        } => Ok(parameters.min_stratigraphic_thickness),
        LayoutDof::Difference {
            kind:
                DifferenceKind::InequalityBelowUpperInterface
                | DifferenceKind::InequalityAboveLowerInterface,
            ..
        } => Ok(0.0),
        _ => Err(AssemblyError::Surfe(Error::InequalityVectorFailure)),
    }
}

fn level(reference: LayoutPointRef, constraints: &Constraints) -> Result<f64, AssemblyError> {
    match reference {
        LayoutPointRef::Interface(index) => constraints
            .interfaces
            .get(index)
            .map(|value| value.level())
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        LayoutPointRef::Inequality(index) => constraints
            .inequalities
            .get(index)
            .map(|value| value.level())
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
    }
}

fn bounded(
    layout: &ConstraintLayout,
    constraints: &Constraints,
    parameters: &Parameters,
    interpolation: &DenseMatrix,
) -> Result<AssemblyConstraints, AssemblyError> {
    let mut bounded_constraints = constraints.clone();
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
            LayoutDof::Difference {
                kind: DifferenceKind::SequencedInterfaces,
                ..
            } => (
                parameters.min_stratigraphic_thickness,
                distance - parameters.min_stratigraphic_thickness,
            ),
            LayoutDof::Difference {
                kind:
                    DifferenceKind::InequalityBelowUpperInterface
                    | DifferenceKind::InequalityAboveLowerInterface,
                ..
            } => (0.0, parameters.min_stratigraphic_thickness),
            LayoutDof::Difference {
                kind: DifferenceKind::SameLevelInterface,
                ..
            } => (
                -parameters.interface_uncertainty,
                2.0 * parameters.interface_uncertainty,
            ),
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
