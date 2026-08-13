//! Frozen Vector Field constraint layout.
//!
//! Source: `surfe_lib/vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`process_input_data`, `get_method_parameters`, `get_interpolation_matrix`,
//! and `get_equality_values`).

use crate::{
    layout::{
        append_section, planar_dofs, ConstraintLayout, IndexRange, LayoutPartitions, LayoutRole,
        LayoutSectionKind, SourceConstraintCounts,
    },
    Constraints, InternalParameters, ModelType, SolverType,
};

pub(crate) fn build(constraints: &Constraints) -> Result<ConstraintLayout, crate::Error> {
    let source_counts = SourceConstraintCounts::from_constraints(constraints);
    let n_planar = constraints.planars.len();
    let n_constraints = 3 * n_planar;
    let internal = InternalParameters {
        n_interface: 0,
        n_planar,
        n_inequality: 0,
        n_tangent: 0,
        n_constraints,
        n_equality: n_constraints,
        modified_basis: false,
        poly_term: false,
        n_poly_terms: 0,
        problem_type: SolverType::Linear,
        restricted_range: false,
    };
    let mut dofs = Vec::with_capacity(n_constraints);
    let mut roles = Vec::with_capacity(n_constraints);
    let mut sections = Vec::new();
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::PlanarDerivatives,
        planar_dofs(n_planar),
        LayoutRole::Equality,
    );
    Ok(ConstraintLayout::new(
        ModelType::VectorField,
        source_counts,
        internal,
        dofs,
        roles,
        sections,
        LayoutPartitions::new(
            IndexRange::new(0, n_constraints),
            IndexRange::new(0, 0),
            IndexRange::new(0, 0),
            IndexRange::new(n_constraints, n_constraints),
        ),
    ))
}
