//! Frozen Lajaunie constraint layout.
//!
//! Source: `surfe_lib/lajaunie.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`_get_increment_pairs`, `process_input_data`, `get_method_parameters`,
//! `get_interpolation_matrix`, `get_equality_values`, and bounded values).

use crate::{
    layout::{
        append_section, declared_polynomial_term_count, planar_dofs, polynomial_dofs, tangent_dofs,
        ConstraintLayout, DifferenceKind, IndexRange, LayoutDof, LayoutPartitions, LayoutPointRef,
        LayoutRole, LayoutSectionKind, SourceConstraintCounts,
    },
    Constraints, Error, InternalParameters, ModelType, Parameters, SolverType,
};

pub(crate) fn build(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<ConstraintLayout, Error> {
    let grouping = constraints
        .interface_grouping()
        .ok_or(Error::NoInterfaceData)?;
    let differences = grouping
        .multi_point_groups()
        .iter()
        .flat_map(|group| {
            let reference = group[0];
            group[1..].iter().map(move |index| LayoutDof::Difference {
                kind: DifferenceKind::SameLevelInterface,
                positive: LayoutPointRef::Interface(reference),
                negative: LayoutPointRef::Interface(*index),
            })
        })
        .collect::<Vec<_>>();
    if differences.is_empty() {
        return Err(Error::NoInterfaceIncrementPairs);
    }

    let source_counts = SourceConstraintCounts::from_constraints(constraints);
    let n_interface = constraints.interfaces.len();
    let n_planar = constraints.planars.len();
    let n_tangent = constraints.tangents.len();
    let n_constraints = differences.len() + 3 * n_planar + n_tangent;
    let restricted = parameters.use_restricted_range;
    let declared_polynomial_terms = declared_polynomial_term_count(parameters.polynomial_order)?
        .checked_sub(1)
        .ok_or(Error::InterpolationMatrixFailure)?;
    let polynomial_terms = if restricted {
        0
    } else {
        declared_polynomial_terms
    };
    let internal = InternalParameters {
        n_interface,
        n_planar,
        n_inequality: 0,
        n_tangent,
        n_constraints,
        n_equality: if restricted { 0 } else { n_constraints },
        modified_basis: restricted,
        poly_term: !restricted,
        n_poly_terms: declared_polynomial_terms,
        problem_type: if restricted {
            SolverType::Quadratic
        } else {
            SolverType::Linear
        },
        restricted_range: restricted,
    };

    let role = if restricted {
        LayoutRole::Bounded
    } else {
        LayoutRole::Equality
    };
    let mut dofs = Vec::with_capacity(n_constraints + polynomial_terms);
    let mut roles = Vec::with_capacity(n_constraints + polynomial_terms);
    let mut sections = Vec::new();
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::SameLevelDifferences,
        differences,
        role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::PlanarDerivatives,
        planar_dofs(n_planar),
        role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::Tangents,
        tangent_dofs(n_tangent),
        role,
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
            IndexRange::new(0, n_constraints),
            IndexRange::new(0, 0),
            IndexRange::new(0, 0),
            IndexRange::new(n_constraints, matrix_size),
        )
    };
    Ok(ConstraintLayout::new(
        ModelType::LajaunieApproach,
        source_counts,
        internal,
        dofs,
        roles,
        sections,
        partitions,
    ))
}
