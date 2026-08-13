//! Frozen Stratigraphic Horizons constraint layout.
//!
//! Source:
//! `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`_get_increment_pairs`, lithostratigraphic pair selection,
//! `process_input_data`, `get_method_parameters`, `get_interpolation_matrix`,
//! `get_equality_values`, and both inequality paths).

use crate::{
    layout::{
        append_section, declared_polynomial_term_count, planar_dofs, tangent_dofs,
        ConstraintLayout, DifferenceKind, IndexRange, LayoutDof, LayoutPartitions, LayoutPointRef,
        LayoutRole, LayoutSectionKind, SourceConstraintCounts,
    },
    Constraints, Error, InterfaceGrouping, InternalParameters, ModelType, Parameters, SolverType,
};

pub(crate) fn build(
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<ConstraintLayout, Error> {
    let grouping = constraints
        .interface_grouping()
        .ok_or(Error::NoInterfaceData)?;
    if constraints.inequalities.iter().any(|inequality| {
        grouping
            .levels_descending()
            .iter()
            .any(|level| *level == inequality.level())
    }) {
        return Err(Error::InvalidInputData);
    }

    let sequenced_interfaces = sequenced_interface_dofs(&grouping);
    let sequenced_inequalities = sequenced_inequality_dofs(constraints, &grouping);
    let same_level = same_level_dofs(&grouping);
    let n_standard_inequality = sequenced_interfaces.len() + sequenced_inequalities.len();
    let n_planar = constraints.planars.len();
    let n_tangent = constraints.tangents.len();
    let n_constraints = n_standard_inequality + same_level.len() + 3 * n_planar + n_tangent;
    let restricted = parameters.use_restricted_range;
    let declared_polynomial_terms = declared_polynomial_term_count(parameters.polynomial_order)?;
    let n_equality = if restricted {
        0
    } else {
        same_level.len() + 3 * n_planar + n_tangent
    };
    let internal = InternalParameters {
        n_interface: constraints.interfaces.len(),
        n_planar,
        n_inequality: if restricted {
            constraints.inequalities.len()
        } else {
            n_standard_inequality
        },
        n_tangent,
        n_constraints,
        n_equality,
        modified_basis: true,
        poly_term: false,
        n_poly_terms: declared_polynomial_terms,
        problem_type: SolverType::Quadratic,
        restricted_range: restricted,
    };
    let source_counts = SourceConstraintCounts::from_constraints(constraints);
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

    let mut dofs = Vec::with_capacity(n_constraints);
    let mut roles = Vec::with_capacity(n_constraints);
    let mut sections = Vec::new();
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::SequencedInterfaceDifferences,
        sequenced_interfaces,
        inequality_role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::SequencedInequalityDifferences,
        sequenced_inequalities,
        inequality_role,
    );
    append_section(
        &mut dofs,
        &mut roles,
        &mut sections,
        LayoutSectionKind::SameLevelDifferences,
        same_level,
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

    let partitions = if restricted {
        LayoutPartitions::new(
            IndexRange::new(0, 0),
            IndexRange::new(0, 0),
            IndexRange::new(0, n_constraints),
            IndexRange::new(n_constraints, n_constraints),
        )
    } else {
        LayoutPartitions::new(
            IndexRange::new(n_standard_inequality, n_constraints),
            IndexRange::new(0, n_standard_inequality),
            IndexRange::new(0, 0),
            IndexRange::new(n_constraints, n_constraints),
        )
    };
    Ok(ConstraintLayout::new(
        ModelType::StratigraphicHorizons,
        source_counts,
        internal,
        dofs,
        roles,
        sections,
        partitions,
    ))
}

fn sequenced_interface_dofs(grouping: &InterfaceGrouping) -> Vec<LayoutDof> {
    grouping
        .reference_indices()
        .windows(2)
        .map(|pair| LayoutDof::Difference {
            kind: DifferenceKind::SequencedInterfaces,
            positive: LayoutPointRef::Interface(pair[0]),
            negative: LayoutPointRef::Interface(pair[1]),
        })
        .collect()
}

fn sequenced_inequality_dofs(
    constraints: &Constraints,
    grouping: &InterfaceGrouping,
) -> Vec<LayoutDof> {
    let levels = grouping.levels_descending();
    let references = grouping.reference_indices();
    let mut dofs = Vec::new();
    for (inequality_index, inequality) in constraints.inequalities.iter().enumerate() {
        let level = inequality.level();
        let above = levels
            .iter()
            .zip(references)
            .filter(|(candidate, _)| **candidate > level)
            .min_by(|(left, _), (right, _)| {
                (**left - level)
                    .partial_cmp(&(**right - level))
                    .expect("constraint levels are finite")
            });
        if let Some((_, reference)) = above {
            dofs.push(LayoutDof::Difference {
                kind: DifferenceKind::InequalityBelowUpperInterface,
                positive: LayoutPointRef::Interface(*reference),
                negative: LayoutPointRef::Inequality(inequality_index),
            });
        }
        let below = levels
            .iter()
            .zip(references)
            .filter(|(candidate, _)| **candidate < level)
            .min_by(|(left, _), (right, _)| {
                (level - **left)
                    .partial_cmp(&(level - **right))
                    .expect("constraint levels are finite")
            });
        if let Some((_, reference)) = below {
            dofs.push(LayoutDof::Difference {
                kind: DifferenceKind::InequalityAboveLowerInterface,
                positive: LayoutPointRef::Inequality(inequality_index),
                negative: LayoutPointRef::Interface(*reference),
            });
        }
    }
    dofs
}

fn same_level_dofs(grouping: &InterfaceGrouping) -> Vec<LayoutDof> {
    grouping
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
        .collect()
}
