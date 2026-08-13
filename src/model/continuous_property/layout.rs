//! Frozen Continuous Property constraint layout.
//!
//! Source:
//! `surfe_lib/continuous_property.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`process_input_data`, `get_method_parameters`, `get_interpolation_matrix`,
//! and `get_equality_values`).

use crate::{
    layout::{
        append_section, ConstraintLayout, IndexRange, LayoutDof, LayoutPartitions, LayoutRole,
        LayoutSectionKind, SourceConstraintCounts,
    },
    Constraints, Error, InternalParameters, ModelType, SolverType,
};

pub(crate) fn build(constraints: &Constraints) -> Result<ConstraintLayout, Error> {
    if constraints.interfaces.is_empty() {
        return Err(Error::NoInterfaceData);
    }
    let source_counts = SourceConstraintCounts::from_constraints(constraints);
    let n_interface = constraints.interfaces.len();
    let internal = InternalParameters {
        n_interface,
        n_planar: 0,
        n_inequality: 0,
        n_tangent: 0,
        n_constraints: n_interface,
        n_equality: n_interface,
        modified_basis: false,
        poly_term: false,
        n_poly_terms: 0,
        problem_type: SolverType::Linear,
        restricted_range: false,
    };
    let mut dofs = Vec::with_capacity(n_interface);
    let mut roles = Vec::with_capacity(n_interface);
    let mut sections = Vec::new();
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::InterfaceValues,
        (0..n_interface).map(|index| LayoutDof::InterfaceValue { index }),
        LayoutRole::Equality,
    );
    Ok(ConstraintLayout::new(
        ModelType::ContinuousProperty,
        source_counts,
        internal,
        dofs,
        roles,
        sections,
        LayoutPartitions::new(
            IndexRange::new(0, n_interface),
            IndexRange::new(0, 0),
            IndexRange::new(0, 0),
            IndexRange::new(n_interface, n_interface),
        ),
    ))
}
