//! Frozen Lajaunie right-hand sides and restricted-range bounds.
//!
//! Source: `surfe_lib/lajaunie.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`get_equality_values` and bounded `get_inequality_values`).

use crate::{
    assembly::normal_component, AssemblyConstraints, AssemblyError, BoundedConstraintSystem,
    ConstraintLayout, Constraints, DenseMatrix, DenseVector, Error, LayoutDof, LayoutPointRef,
    Parameters,
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
        Ok(AssemblyConstraints::Linear {
            right_hand_side: DenseVector::from_values(
                layout
                    .dofs()
                    .iter()
                    .map(|dof| equality_value(dof, constraints))
                    .collect::<Result<Vec<_>, _>>()?,
            ),
        })
    }
}

fn equality_value(dof: &LayoutDof, constraints: &Constraints) -> Result<f64, AssemblyError> {
    match dof {
        LayoutDof::Difference {
            positive, negative, ..
        } => Ok(level(*positive, constraints)? - level(*negative, constraints)?),
        LayoutDof::PlanarDerivative { index, axis } => constraints
            .planars
            .get(*index)
            .map(|value| normal_component(value.normal(), *axis))
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        LayoutDof::Tangent { .. } | LayoutDof::PolynomialTerm { .. } => Ok(0.0),
        _ => Err(AssemblyError::Surfe(Error::EqualityVectorFailure)),
    }
}

fn level(reference: LayoutPointRef, constraints: &Constraints) -> Result<f64, AssemblyError> {
    match reference {
        LayoutPointRef::Interface(index) => constraints
            .interfaces
            .get(index)
            .map(|value| value.level())
            .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        LayoutPointRef::Inequality(_) => Err(AssemblyError::Surfe(Error::EqualityVectorFailure)),
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
    let mut lower = Vec::with_capacity(layout.matrix_size());
    let mut range = Vec::with_capacity(layout.matrix_size());
    for dof in layout.dofs() {
        let bounds = match dof {
            LayoutDof::Difference { .. } => (
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
