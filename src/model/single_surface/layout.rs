//! Frozen Single Surface constraint layout.
//!
//! Source: `surfe_lib/single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`process_input_data`, `get_method_parameters`, `get_interpolation_matrix`,
//! `get_equality_values`, and both inequality paths).

use crate::{
    layout::{
        append_section, declared_polynomial_term_count, planar_dofs, polynomial_dofs, tangent_dofs,
        ConstraintLayout, IndexRange, LayoutDof, LayoutPartitions, LayoutRole, LayoutSectionKind,
        SourceConstraintCounts,
    },
    Constraints, Error, InternalParameters, ModelType, Parameters, SolverType,
};

pub(crate) fn build(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<ConstraintLayout, Error> {
    if constraints.interfaces.is_empty() {
        return Err(Error::NoInterfaceData);
    }

    let source_counts = SourceConstraintCounts::from_constraints(constraints);
    let n_inequality = constraints.inequalities.len();
    let n_interface = constraints.interfaces.len();
    let n_planar = constraints.planars.len();
    let n_tangent = constraints.tangents.len();
    let n_constraints = n_inequality + n_interface + 3 * n_planar + n_tangent;
    let restricted = parameters.use_restricted_range;
    let quadratic = n_inequality != 0 || restricted;
    let n_equality = if restricted {
        0
    } else {
        n_interface + 3 * n_planar + n_tangent
    };
    let declared_polynomial_terms = declared_polynomial_term_count(parameters.polynomial_order)?;
    let polynomial_terms = if quadratic {
        0
    } else {
        declared_polynomial_terms
    };
    let internal = InternalParameters {
        n_interface,
        n_planar,
        n_inequality,
        n_tangent,
        n_constraints,
        n_equality,
        modified_basis: quadratic,
        poly_term: !quadratic,
        n_poly_terms: declared_polynomial_terms,
        problem_type: if quadratic {
            SolverType::Quadratic
        } else {
            SolverType::Linear
        },
        restricted_range: restricted,
    };

    let mut dofs = Vec::with_capacity(n_constraints + polynomial_terms);
    let mut roles = Vec::with_capacity(n_constraints + polynomial_terms);
    let mut sections = Vec::new();
    let inequality_role = if restricted {
        LayoutRole::Bounded
    } else {
        LayoutRole::Inequality
    };
    let equality_role = if restricted {
        LayoutRole::Bounded
    } else {
        LayoutRole::Equality
    };
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::InequalityValues,
        (0..n_inequality).map(|index| LayoutDof::InequalityValue { index }),
        inequality_role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::InterfaceValues,
        (0..n_interface).map(|index| LayoutDof::InterfaceValue { index }),
        equality_role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::PlanarDerivatives,
        planar_dofs(n_planar),
        equality_role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::Tangents,
        tangent_dofs(n_tangent),
        equality_role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::Polynomial,
        polynomial_dofs(polynomial_terms),
        LayoutRole::Polynomial,
    );

    let matrix_size = n_constraints + polynomial_terms;
    let partitions = if restricted {
        LayoutPartitions::new(
            IndexRange::new(0, 0),
            IndexRange::new(0, 0),
            IndexRange::new(0, n_constraints),
            IndexRange::new(n_constraints, matrix_size),
        )
    } else {
        LayoutPartitions::new(
            IndexRange::new(n_inequality, n_constraints),
            IndexRange::new(0, n_inequality),
            IndexRange::new(0, 0),
            IndexRange::new(n_constraints, matrix_size),
        )
    };
    Ok(ConstraintLayout::new(
        ModelType::SingleSurface,
        source_counts,
        internal,
        dofs,
        roles,
        sections,
        partitions,
    ))
}
