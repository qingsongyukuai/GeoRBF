use georbf::{
    constraint_layout, Axis, ConstraintLayout, Constraints, DifferenceKind, Error, IndexRange,
    Inequality, Interface, InternalParameters, LayoutDof, LayoutPointRef, LayoutRole,
    LayoutSectionKind, ModelType, Parameters, Planar, PolynomialOrder, SolverType, Tangent,
};

fn interface(x: f64, level: f64) -> Interface {
    Interface::new(x, 0.0, 0.0, level).unwrap()
}

fn inequality(x: f64, level: f64) -> Inequality {
    Inequality::new(x, 1.0, 0.0, level).unwrap()
}

fn planar(x: f64) -> Planar {
    Planar::from_normal(x, 2.0, 0.0, 0.0, 0.0, 1.0).unwrap()
}

fn tangent(x: f64) -> Tangent {
    Tangent::new(x, 3.0, 0.0, 1.0, 0.0, 0.0).unwrap()
}

fn representative_constraints() -> Constraints {
    Constraints {
        inequalities: vec![inequality(20.0, 25.0), inequality(21.0, 5.0)],
        interfaces: vec![
            interface(0.0, 30.0),
            interface(1.0, 30.0),
            interface(2.0, 20.0),
            interface(3.0, 10.0),
            interface(4.0, 10.0),
        ],
        planars: vec![planar(30.0), planar(31.0)],
        tangents: vec![tangent(40.0)],
    }
}

fn parameters(model: ModelType) -> Parameters {
    Parameters {
        model_type: model,
        polynomial_order: 1,
        ..Parameters::default()
    }
}

fn layout(model: ModelType) -> ConstraintLayout {
    constraint_layout(model, &representative_constraints(), &parameters(model)).unwrap()
}

#[test]
fn five_models_match_frozen_dimensions_sections_and_partitions() {
    let single = layout(ModelType::SingleSurface);
    assert_eq!(single.matrix_size(), 14);
    assert_eq!(single.constraint_dof_count(), 14);
    assert_eq!(single.polynomial_dof_count(), 0);
    assert_eq!(
        single.internal_parameters(),
        &InternalParameters {
            n_interface: 5,
            n_planar: 2,
            n_inequality: 2,
            n_tangent: 1,
            n_constraints: 14,
            n_equality: 12,
            modified_basis: true,
            poly_term: false,
            n_poly_terms: 4,
            problem_type: SolverType::Quadratic,
            restricted_range: false,
        }
    );
    assert_eq!(single.partitions().inequality(), IndexRange::new(0, 2));
    assert_eq!(single.partitions().equality(), IndexRange::new(2, 14));
    assert_eq!(single.partitions().bounded(), IndexRange::new(0, 0));
    assert_eq!(single.partitions().polynomial(), IndexRange::new(14, 14));
    assert_eq!(
        single.section(LayoutSectionKind::InequalityValues),
        Some(IndexRange::new(0, 2))
    );
    assert_eq!(
        single.section(LayoutSectionKind::InterfaceValues),
        Some(IndexRange::new(2, 7))
    );
    assert_eq!(
        single.section(LayoutSectionKind::PlanarDerivatives),
        Some(IndexRange::new(7, 13))
    );
    assert_eq!(
        single.section(LayoutSectionKind::Tangents),
        Some(IndexRange::new(13, 14))
    );

    let lajaunie = layout(ModelType::LajaunieApproach);
    assert_eq!(lajaunie.matrix_size(), 12);
    assert_eq!(lajaunie.constraint_dof_count(), 9);
    assert_eq!(lajaunie.polynomial_dof_count(), 3);
    assert_eq!(
        lajaunie.internal_parameters(),
        &InternalParameters {
            n_interface: 5,
            n_planar: 2,
            n_inequality: 0,
            n_tangent: 1,
            n_constraints: 9,
            n_equality: 9,
            modified_basis: false,
            poly_term: true,
            n_poly_terms: 3,
            problem_type: SolverType::Linear,
            restricted_range: false,
        }
    );
    assert_eq!(lajaunie.partitions().equality(), IndexRange::new(0, 9));
    assert_eq!(lajaunie.partitions().polynomial(), IndexRange::new(9, 12));
    assert_eq!(
        lajaunie.section(LayoutSectionKind::SameLevelDifferences),
        Some(IndexRange::new(0, 2))
    );
    assert_eq!(
        lajaunie.section(LayoutSectionKind::PlanarDerivatives),
        Some(IndexRange::new(2, 8))
    );

    let stratigraphic = layout(ModelType::StratigraphicHorizons);
    assert_eq!(stratigraphic.matrix_size(), 14);
    assert_eq!(stratigraphic.constraint_dof_count(), 14);
    assert_eq!(stratigraphic.polynomial_dof_count(), 0);
    assert_eq!(
        stratigraphic.internal_parameters(),
        &InternalParameters {
            n_interface: 5,
            n_planar: 2,
            n_inequality: 5,
            n_tangent: 1,
            n_constraints: 14,
            n_equality: 9,
            modified_basis: true,
            poly_term: false,
            n_poly_terms: 4,
            problem_type: SolverType::Quadratic,
            restricted_range: false,
        }
    );
    assert_eq!(
        stratigraphic.partitions().inequality(),
        IndexRange::new(0, 5)
    );
    assert_eq!(
        stratigraphic.partitions().equality(),
        IndexRange::new(5, 14)
    );
    assert_eq!(
        stratigraphic.section(LayoutSectionKind::SequencedInterfaceDifferences),
        Some(IndexRange::new(0, 2))
    );
    assert_eq!(
        stratigraphic.section(LayoutSectionKind::SequencedInequalityDifferences),
        Some(IndexRange::new(2, 5))
    );
    assert_eq!(
        stratigraphic.section(LayoutSectionKind::SameLevelDifferences),
        Some(IndexRange::new(5, 7))
    );

    let continuous = layout(ModelType::ContinuousProperty);
    assert_eq!(continuous.matrix_size(), 5);
    assert_eq!(
        continuous.internal_parameters(),
        &InternalParameters {
            n_interface: 5,
            n_planar: 0,
            n_inequality: 0,
            n_tangent: 0,
            n_constraints: 5,
            n_equality: 5,
            modified_basis: false,
            poly_term: false,
            n_poly_terms: 0,
            problem_type: SolverType::Linear,
            restricted_range: false,
        }
    );
    assert_eq!(continuous.partitions().equality(), IndexRange::new(0, 5));
    assert_eq!(
        continuous.section(LayoutSectionKind::InterfaceValues),
        Some(IndexRange::new(0, 5))
    );
    assert_eq!(continuous.internal_parameters().n_planar, 0);
    assert_eq!(continuous.internal_parameters().n_tangent, 0);

    let vector = layout(ModelType::VectorField);
    assert_eq!(vector.matrix_size(), 6);
    assert_eq!(
        vector.internal_parameters(),
        &InternalParameters {
            n_interface: 0,
            n_planar: 2,
            n_inequality: 0,
            n_tangent: 0,
            n_constraints: 6,
            n_equality: 6,
            modified_basis: false,
            poly_term: false,
            n_poly_terms: 0,
            problem_type: SolverType::Linear,
            restricted_range: false,
        }
    );
    assert_eq!(vector.partitions().equality(), IndexRange::new(0, 6));
    assert_eq!(
        vector.section(LayoutSectionKind::PlanarDerivatives),
        Some(IndexRange::new(0, 6))
    );
    assert_eq!(vector.internal_parameters().n_interface, 0);
}

#[test]
fn row_and_column_dof_labels_preserve_each_models_source_order() {
    let single = layout(ModelType::SingleSurface);
    assert_eq!(
        single.dof(0),
        Some(&LayoutDof::InequalityValue { index: 0 })
    );
    assert_eq!(single.dof(2), Some(&LayoutDof::InterfaceValue { index: 0 }));
    assert_eq!(
        single.dof(7),
        Some(&LayoutDof::PlanarDerivative {
            index: 0,
            axis: Axis::X,
        })
    );
    assert_eq!(
        single.dof(10),
        Some(&LayoutDof::PlanarDerivative {
            index: 1,
            axis: Axis::X,
        })
    );
    assert_eq!(single.dof(13), Some(&LayoutDof::Tangent { index: 0 }));

    let lajaunie = layout(ModelType::LajaunieApproach);
    let same_level_30 = LayoutDof::Difference {
        kind: DifferenceKind::SameLevelInterface,
        positive: LayoutPointRef::Interface(0),
        negative: LayoutPointRef::Interface(1),
    };
    let same_level_10 = LayoutDof::Difference {
        kind: DifferenceKind::SameLevelInterface,
        positive: LayoutPointRef::Interface(3),
        negative: LayoutPointRef::Interface(4),
    };
    assert_eq!(lajaunie.dof(0), Some(&same_level_30));
    assert_eq!(lajaunie.dof(1), Some(&same_level_10));
    assert_eq!(lajaunie.index_of(&same_level_10), Some(1));
    assert_eq!(lajaunie.role(9), Some(LayoutRole::Polynomial));
    assert_eq!(
        lajaunie.dof(9),
        Some(&LayoutDof::PolynomialTerm { index: 0 })
    );

    let stratigraphic = layout(ModelType::StratigraphicHorizons);
    assert_eq!(
        stratigraphic.dof(0),
        Some(&LayoutDof::Difference {
            kind: DifferenceKind::SequencedInterfaces,
            positive: LayoutPointRef::Interface(0),
            negative: LayoutPointRef::Interface(2),
        })
    );
    assert_eq!(
        stratigraphic.dof(2),
        Some(&LayoutDof::Difference {
            kind: DifferenceKind::InequalityBelowUpperInterface,
            positive: LayoutPointRef::Interface(0),
            negative: LayoutPointRef::Inequality(0),
        })
    );
    assert_eq!(
        stratigraphic.dof(3),
        Some(&LayoutDof::Difference {
            kind: DifferenceKind::InequalityAboveLowerInterface,
            positive: LayoutPointRef::Inequality(0),
            negative: LayoutPointRef::Interface(2),
        })
    );
    assert_eq!(
        stratigraphic.dof(4),
        Some(&LayoutDof::Difference {
            kind: DifferenceKind::InequalityBelowUpperInterface,
            positive: LayoutPointRef::Interface(3),
            negative: LayoutPointRef::Inequality(1),
        })
    );
    assert_eq!(stratigraphic.role(4), Some(LayoutRole::Inequality));
    assert_eq!(stratigraphic.role(5), Some(LayoutRole::Equality));
}

#[test]
fn polynomial_and_restricted_range_layouts_follow_frozen_solver_branches() {
    let mut constraints = representative_constraints();
    constraints.inequalities.clear();

    for (order, expected_terms) in [
        (PolynomialOrder::Zero, 1usize),
        (PolynomialOrder::First, 4),
        (PolynomialOrder::Second, 10),
    ] {
        let mut params = parameters(ModelType::SingleSurface);
        params.polynomial_order = match order {
            PolynomialOrder::Zero => 0,
            PolynomialOrder::First => 1,
            PolynomialOrder::Second => 2,
        };
        let layout = constraint_layout(ModelType::SingleSurface, &constraints, &params).unwrap();
        assert_eq!(layout.constraint_dof_count(), 12);
        assert_eq!(layout.polynomial_dof_count(), expected_terms);
        assert_eq!(layout.matrix_size(), 12 + expected_terms);
        assert_eq!(layout.partitions().equality(), IndexRange::new(0, 12));
        assert_eq!(
            layout.partitions().polynomial(),
            IndexRange::new(12, 12 + expected_terms)
        );
    }

    // Frozen `get_method_parameters` computes the declared tetrahedral count
    // before later polynomial-basis validation. T16 records that layout only.
    let mut unsupported_order = parameters(ModelType::SingleSurface);
    unsupported_order.polynomial_order = 3;
    let unsupported =
        constraint_layout(ModelType::SingleSurface, &constraints, &unsupported_order).unwrap();
    assert_eq!(unsupported.polynomial_dof_count(), 20);
    assert_eq!(unsupported.matrix_size(), 32);

    let mut lajaunie_params = parameters(ModelType::LajaunieApproach);
    lajaunie_params.polynomial_order = 2;
    let lajaunie =
        constraint_layout(ModelType::LajaunieApproach, &constraints, &lajaunie_params).unwrap();
    assert_eq!(lajaunie.polynomial_dof_count(), 9);

    for model in [
        ModelType::SingleSurface,
        ModelType::LajaunieApproach,
        ModelType::StratigraphicHorizons,
    ] {
        let mut params = parameters(model);
        params.use_restricted_range = true;
        let layout = constraint_layout(model, &representative_constraints(), &params).unwrap();
        assert_eq!(layout.polynomial_dof_count(), 0);
        assert_eq!(
            layout.partitions().bounded(),
            IndexRange::new(0, layout.matrix_size())
        );
        assert_eq!(layout.partitions().equality(), IndexRange::new(0, 0));
        assert_eq!(layout.partitions().inequality(), IndexRange::new(0, 0));
        assert!(layout.internal_parameters().modified_basis);
        assert!(layout.internal_parameters().restricted_range);
    }
}

#[test]
fn empty_categories_and_invalid_groupings_match_frozen_errors() {
    let empty = Constraints::default();
    for model in [
        ModelType::SingleSurface,
        ModelType::LajaunieApproach,
        ModelType::StratigraphicHorizons,
        ModelType::ContinuousProperty,
    ] {
        assert_eq!(
            constraint_layout(model, &empty, &parameters(model)),
            Err(Error::NoInterfaceData)
        );
    }

    let vector = constraint_layout(
        ModelType::VectorField,
        &empty,
        &parameters(ModelType::VectorField),
    )
    .unwrap();
    assert_eq!(vector.matrix_size(), 0);
    assert_eq!(
        vector.section(LayoutSectionKind::PlanarDerivatives),
        Some(IndexRange::new(0, 0))
    );

    let singleton = Constraints {
        interfaces: vec![interface(0.0, 10.0), interface(1.0, 20.0)],
        ..Constraints::default()
    };
    assert_eq!(
        constraint_layout(
            ModelType::LajaunieApproach,
            &singleton,
            &parameters(ModelType::LajaunieApproach)
        ),
        Err(Error::NoInterfaceIncrementPairs)
    );
    let stratigraphic = constraint_layout(
        ModelType::StratigraphicHorizons,
        &singleton,
        &parameters(ModelType::StratigraphicHorizons),
    )
    .unwrap();
    assert_eq!(stratigraphic.matrix_size(), 1);

    let invalid_stratigraphic = Constraints {
        inequalities: vec![inequality(4.0, 20.0)],
        ..singleton
    };
    assert_eq!(
        constraint_layout(
            ModelType::StratigraphicHorizons,
            &invalid_stratigraphic,
            &parameters(ModelType::StratigraphicHorizons)
        ),
        Err(Error::InvalidInputData)
    );
}

#[test]
fn cleaned_input_permutations_produce_identical_layouts_for_all_models() {
    let original = representative_constraints();
    let mut forward = original.clone();
    forward.interfaces.push(original.interfaces[0].clone());
    forward.planars.push(original.planars[0].clone());
    forward.inequalities.push(original.inequalities[0].clone());
    forward.tangents.push(original.tangents[0].clone());

    let mut reversed = forward.clone();
    reversed.interfaces.reverse();
    reversed.planars.reverse();
    reversed.inequalities.reverse();
    reversed.tangents.reverse();

    forward.remove_collocated();
    reversed.remove_collocated();

    for model in [
        ModelType::SingleSurface,
        ModelType::LajaunieApproach,
        ModelType::StratigraphicHorizons,
        ModelType::ContinuousProperty,
        ModelType::VectorField,
    ] {
        let forward_layout = constraint_layout(model, &forward, &parameters(model)).unwrap();
        let reverse_layout = constraint_layout(model, &reversed, &parameters(model)).unwrap();
        assert_eq!(forward_layout, reverse_layout, "{model:?}");
    }
}

#[test]
fn debug_snapshot_is_exact_stable_and_self_describing() {
    let constraints = Constraints {
        planars: vec![planar(1.0), planar(2.0)],
        ..Constraints::default()
    };
    let layout = constraint_layout(
        ModelType::VectorField,
        &constraints,
        &parameters(ModelType::VectorField),
    )
    .unwrap();
    assert_eq!(
        layout.debug_snapshot(),
        "model=Vector_field solver=Linear modified=false restricted=false\n\
source inequality=0 interface=0 planar=2 tangent=0\n\
internal inequality=0 interface=0 planar=2 tangent=0 constraints=6 equality=6 polynomial=0\n\
matrix size=6 equality=0..6 inequality=0..0 bounded=0..0 polynomial=6..6\n\
section planar_derivatives=0..6\n\
0 equality planar[0].dx\n\
1 equality planar[0].dy\n\
2 equality planar[0].dz\n\
3 equality planar[1].dx\n\
4 equality planar[1].dy\n\
5 equality planar[1].dz"
    );
}
