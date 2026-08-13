//! Frozen Vector Field normal-component right-hand side.
//!
//! Source: `surfe_lib/vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`get_equality_values`).

use crate::{
    assembly::normal_component, AssemblyConstraints, AssemblyError, ConstraintLayout, Constraints,
    DenseVector, Error, LayoutDof,
};

pub(crate) fn build(
    layout: &ConstraintLayout,
    constraints: &Constraints,
) -> Result<AssemblyConstraints, AssemblyError> {
    let values = layout
        .dofs()
        .iter()
        .map(|dof| match dof {
            LayoutDof::PlanarDerivative { index, axis } => constraints
                .planars
                .get(*index)
                .map(|value| normal_component(value.normal(), *axis))
                .ok_or(AssemblyError::Surfe(Error::EqualityVectorFailure)),
            _ => Err(AssemblyError::Surfe(Error::EqualityVectorFailure)),
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(AssemblyConstraints::Linear {
        right_hand_side: DenseVector::from_values(values),
    })
}
