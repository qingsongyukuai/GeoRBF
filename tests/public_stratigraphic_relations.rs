use georbf::fit::BoundActiveState;
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::kernel::FieldEnergyNormalization;
use georbf::observation::{FieldValueObservation, GradientObservation, QuadraticPenalty};
use georbf::problem::{BuildError, BuilderConfigurationError};
use georbf::relation::{
    AdditiveFieldGauge, FieldLevelOrder, FieldValueBound, HorizonBuilder, LinearViolationPenalty,
    MinimumFieldSeparation, MinimumFieldSeparationError, OlderThan, SharedLevelSetRelationKind,
    StratigraphicFieldDirection, YoungerThan,
};
use georbf::{GroupId, Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).unwrap()
}

fn builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["x", "y", "z"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        FieldUnitLabel::new("stratigraphic-unit"),
    )
}

fn horizon(group: &str, source: &str, z: f64) -> georbf::relation::Horizon {
    let mut horizon = HorizonBuilder::new(GroupId::new(group));
    horizon
        .add_member(SourceId::new(source), point(0.0, 0.0, z))
        .unwrap();
    horizon.build().unwrap()
}

#[test]
fn minimum_field_separation_is_a_checked_positive_field_quantity() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            MinimumFieldSeparation::try_new(value),
            Err(MinimumFieldSeparationError::NotFinite)
        );
    }
    for value in [-1.0, -0.0, 0.0] {
        assert_eq!(
            MinimumFieldSeparation::try_new(value),
            Err(MinimumFieldSeparationError::NotPositive)
        );
    }
    assert_eq!(MinimumFieldSeparation::try_new(2.5).unwrap().value(), 2.5);
}

#[test]
fn age_relations_require_one_explicit_atomic_field_direction_but_order_does_not() {
    let younger = GroupId::new("younger");
    let older = GroupId::new("older");
    let separation = MinimumFieldSeparation::try_new(1.0).unwrap();
    let mut problem = builder();
    problem
        .add(YoungerThan::hard(
            SourceId::new("age"),
            younger.clone(),
            older.clone(),
            separation,
        ))
        .unwrap();
    problem
        .add(FieldLevelOrder::hard(
            SourceId::new("order"),
            older.clone(),
            younger.clone(),
        ))
        .unwrap();
    problem.add(horizon("older", "older/member", 0.0)).unwrap();
    problem
        .add(horizon("younger", "younger/member", 1.0))
        .unwrap();

    let failure = problem.build().unwrap_err();
    assert_eq!(
        failure.errors(),
        &[BuildError::MissingStratigraphicFieldDirection]
    );
    let mut problem = failure.into_builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    assert_eq!(
        problem.set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger),
        Err(BuilderConfigurationError::StratigraphicFieldDirectionAlreadySet)
    );
    assert_eq!(
        problem.set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardOlder),
        Err(BuilderConfigurationError::StratigraphicFieldDirectionAlreadySet)
    );
    let snapshot = problem.build().unwrap();
    assert_eq!(
        snapshot.stratigraphic_field_direction(),
        Some(StratigraphicFieldDirection::TowardYounger)
    );
    assert_eq!(snapshot.stratigraphic_age_relation_count(), 1);
    assert_eq!(snapshot.field_level_order_count(), 1);
}

#[test]
fn shared_level_set_relations_allow_forward_references_and_sort_dangling_groups() {
    let separation = MinimumFieldSeparation::try_new(1.0).unwrap();
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardOlder)
        .unwrap();
    problem
        .add(OlderThan::hard(
            SourceId::new("b-relation"),
            GroupId::new("missing-z"),
            GroupId::new("present"),
            separation,
        ))
        .unwrap();
    problem
        .add(FieldLevelOrder::hard(
            SourceId::new("a-relation"),
            GroupId::new("missing-a"),
            GroupId::new("present"),
        ))
        .unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(
                SourceId::new("z-gauge"),
                GroupId::new("missing-gauge"),
                0.0,
            )
            .unwrap(),
        )
        .unwrap();
    problem
        .add(horizon("present", "present/member", 0.0))
        .unwrap();

    let failure = problem.build().unwrap_err();
    assert_eq!(
        failure.errors(),
        &[
            BuildError::UnknownGroupReference {
                source_id: SourceId::new("a-relation"),
                group_id: GroupId::new("missing-a"),
            },
            BuildError::UnknownGroupReference {
                source_id: SourceId::new("b-relation"),
                group_id: GroupId::new("missing-z"),
            },
            BuildError::UnknownGroupReference {
                source_id: SourceId::new("z-gauge"),
                group_id: GroupId::new("missing-gauge"),
            },
        ]
    );
}

#[test]
fn hard_self_reverse_and_strict_cycles_are_direct_input_conflicts() {
    let separation = MinimumFieldSeparation::try_new(1.0).unwrap();
    let cases = [
        vec![YoungerThan::hard(
            SourceId::new("self"),
            GroupId::new("a"),
            GroupId::new("a"),
            separation,
        )],
        vec![
            YoungerThan::hard(
                SourceId::new("a-before-b"),
                GroupId::new("b"),
                GroupId::new("a"),
                separation,
            ),
            YoungerThan::hard(
                SourceId::new("b-before-a"),
                GroupId::new("a"),
                GroupId::new("b"),
                separation,
            ),
        ],
        vec![
            YoungerThan::hard(
                SourceId::new("a-to-b"),
                GroupId::new("b"),
                GroupId::new("a"),
                separation,
            ),
            YoungerThan::hard(
                SourceId::new("b-to-c"),
                GroupId::new("c"),
                GroupId::new("b"),
                separation,
            ),
            YoungerThan::hard(
                SourceId::new("c-to-a"),
                GroupId::new("a"),
                GroupId::new("c"),
                separation,
            ),
        ],
    ];

    for relations in cases {
        let mut problem = builder();
        problem
            .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
            .unwrap();
        let mut groups = relations
            .iter()
            .flat_map(|relation| {
                [
                    relation.younger_group_id().clone(),
                    relation.older_group_id().clone(),
                ]
            })
            .collect::<Vec<_>>();
        groups.sort();
        groups.dedup();
        for (index, group) in groups.iter().enumerate() {
            problem
                .add(horizon(
                    group.as_str(),
                    &format!("{}/member", group.as_str()),
                    index as f64,
                ))
                .unwrap();
        }
        problem
            .add(
                AdditiveFieldGauge::at_level_set(SourceId::new("gauge"), groups[0].clone(), 0.0)
                    .unwrap(),
            )
            .unwrap();
        for relation in relations {
            problem.add(relation).unwrap();
        }

        let failure = problem.build().unwrap().fit().unwrap_err();
        assert_eq!(
            failure.diagnosis(),
            georbf::diagnostics::ProblemDiagnosis::DirectInputConflict
        );
        let evidence = failure
            .report()
            .shared_level_set_relation_conflict()
            .expect("the graph proof is retained");
        assert!(!evidence.source_ids().is_empty());
        assert_eq!(evidence.source_ids().len(), evidence.semantic_roles().len());
        assert!(!evidence.group_ids().is_empty());
        assert!(!evidence.backend_invoked());
        assert!(failure.report().attempts().is_empty());
        let witness = failure
            .report()
            .conflict_witness()
            .expect("a strict relation cycle has a canonical source witness");
        assert_eq!(witness.source_ids(), evidence.source_ids());
        assert_eq!(witness.canonical_residual(), 0.0);
        assert!(witness.separation_margin() >= witness.separation_limit());
        assert!(witness.provenance_verified());
        assert!(!witness.backend_invoked());
    }
}

#[test]
fn a_strict_chain_conflicting_with_non_strict_field_order_retains_the_full_proof() {
    let separation = MinimumFieldSeparation::try_new(2.0).unwrap();
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    for (group, z) in [("a", 0.0), ("b", 1.0), ("c", 2.0)] {
        problem
            .add(horizon(group, &format!("{group}/member"), z))
            .unwrap();
    }
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("gauge"), GroupId::new("a"), 0.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(YoungerThan::hard(
            SourceId::new("a-to-b"),
            GroupId::new("b"),
            GroupId::new("a"),
            separation,
        ))
        .unwrap();
    problem
        .add(OlderThan::hard(
            SourceId::new("b-to-c"),
            GroupId::new("b"),
            GroupId::new("c"),
            separation,
        ))
        .unwrap();
    problem
        .add(FieldLevelOrder::hard(
            SourceId::new("c-no-greater-than-a"),
            GroupId::new("c"),
            GroupId::new("a"),
        ))
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    assert_eq!(
        failure.diagnosis(),
        georbf::diagnostics::ProblemDiagnosis::DirectInputConflict
    );
    let evidence = failure
        .report()
        .shared_level_set_relation_conflict()
        .unwrap();
    assert_eq!(
        evidence
            .source_ids()
            .iter()
            .map(SourceId::as_str)
            .collect::<Vec<_>>(),
        ["a-to-b", "b-to-c", "c-no-greater-than-a"]
    );
    assert_eq!(
        evidence
            .group_ids()
            .iter()
            .map(GroupId::as_str)
            .collect::<Vec<_>>(),
        ["a", "b", "c"]
    );
    assert_eq!(evidence.semantic_roles().len(), 3);
}

#[test]
fn younger_than_fits_and_recovers_physical_observables_in_both_field_directions() {
    for (direction, gradient_z, expected_younger) in [
        (StratigraphicFieldDirection::TowardYounger, 1.0, 2.0),
        (StratigraphicFieldDirection::TowardOlder, -1.0, -2.0),
    ] {
        let mut problem = builder();
        problem
            .set_stratigraphic_field_direction(direction)
            .unwrap();
        problem.add(horizon("older", "older/member", 0.0)).unwrap();
        problem
            .add(horizon("younger", "younger/member", 2.0))
            .unwrap();
        problem
            .add(
                AdditiveFieldGauge::at_level_set(
                    SourceId::new("gauge"),
                    GroupId::new("older"),
                    0.0,
                )
                .unwrap(),
            )
            .unwrap();
        problem
            .add(GradientObservation::new(
                SourceId::new("gradient"),
                point(0.0, 0.0, 1.0),
                Vector3::try_new(0.0, 0.0, gradient_z).unwrap(),
            ))
            .unwrap();
        problem
            .add(YoungerThan::hard(
                SourceId::new("age"),
                GroupId::new("younger"),
                GroupId::new("older"),
                MinimumFieldSeparation::try_new(1.0).unwrap(),
            ))
            .unwrap();

        let success = problem.build().unwrap().fit().unwrap();
        let younger = success
            .model()
            .shared_level_value(&GroupId::new("younger"))
            .unwrap();
        let older = success
            .model()
            .shared_level_value(&GroupId::new("older"))
            .unwrap();
        assert!((older - 0.0).abs() <= 1.0e-7);
        assert!(
            (younger - expected_younger).abs() <= 1.0e-7,
            "direction={direction:?}, younger={younger:e}, expected={expected_younger:e}, tolerance=1e-7"
        );
        let query = success.model().evaluate(point(0.0, 0.0, 1.0)).unwrap();
        assert!((query.value() - gradient_z).abs() <= 1.0e-7);
        assert!((query.gradient().components()[2] - gradient_z).abs() <= 1.0e-7);

        let assessment = &success.report().shared_level_set_relations()[0];
        assert_eq!(assessment.source_id().as_str(), "age");
        assert_eq!(assessment.kind(), SharedLevelSetRelationKind::YoungerThan);
        assert_eq!(assessment.younger_group_id().unwrap().as_str(), "younger");
        assert_eq!(assessment.older_group_id().unwrap().as_str(), "older");
        assert_eq!(assessment.lower_group_id(), None);
        assert_eq!(assessment.upper_group_id(), None);
        assert!((assessment.recovered_younger_value().unwrap() - younger).abs() <= 1.0e-7);
        assert!((assessment.recovered_older_value().unwrap() - older).abs() <= 1.0e-7);
        assert_eq!(assessment.recovered_lower_value(), None);
        assert_eq!(assessment.recovered_upper_value(), None);
        assert_eq!(assessment.field_direction(), Some(direction));
        assert_eq!(assessment.minimum_separation().unwrap().value(), 1.0);
        assert!((assessment.recovered_field_separation() - 2.0).abs() <= 1.0e-7);
        assert!((assessment.slack() - 1.0).abs() <= 1.0e-7);
        assert_eq!(assessment.violation(), 0.0);
        assert_eq!(assessment.active_state(), BoundActiveState::Inactive);
        assert!(assessment.loss().is_none());
        assert!(!success.report().attempts().is_empty());
    }
}

#[test]
fn field_level_order_is_non_strict_and_allows_equal_shared_values_without_age_direction() {
    let mut problem = builder();
    problem.add(horizon("a", "a/member", 0.0)).unwrap();
    problem.add(horizon("b", "b/member", 0.0)).unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("gauge"), GroupId::new("a"), 3.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.0, 0.0, 0.0),
            Vector3::try_new(1.0, 0.0, 0.0).unwrap(),
        ))
        .unwrap();
    problem
        .add(FieldLevelOrder::hard(
            SourceId::new("order"),
            GroupId::new("a"),
            GroupId::new("b"),
        ))
        .unwrap();

    let success = problem.build().unwrap().fit().unwrap();
    let assessment = &success.report().shared_level_set_relations()[0];
    assert_eq!(
        assessment.kind(),
        SharedLevelSetRelationKind::FieldLevelOrder
    );
    assert_eq!(assessment.field_direction(), None);
    assert_eq!(assessment.minimum_separation(), None);
    assert_eq!(assessment.younger_group_id(), None);
    assert_eq!(assessment.older_group_id(), None);
    assert_eq!(assessment.lower_group_id().unwrap().as_str(), "a");
    assert_eq!(assessment.upper_group_id().unwrap().as_str(), "b");
    assert!((assessment.recovered_lower_value().unwrap() - 3.0).abs() <= 1.0e-8);
    assert!((assessment.recovered_upper_value().unwrap() - 3.0).abs() <= 1.0e-8);
    assert_eq!(assessment.recovered_younger_value(), None);
    assert_eq!(assessment.recovered_older_value(), None);
    assert!(assessment.recovered_field_separation().abs() <= 1.0e-8);
    assert!(assessment.slack() <= assessment.tolerance());
    assert_eq!(assessment.active_state(), BoundActiveState::Active);
}

#[test]
fn soft_age_duplicates_keep_independent_quadratic_and_linear_violations() {
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    problem
        .set_field_energy_normalization(FieldEnergyNormalization::try_new(1.0).unwrap())
        .unwrap();
    problem.add(horizon("older", "older/member", 0.0)).unwrap();
    problem
        .add(horizon("younger", "younger/member", 0.0))
        .unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("gauge"), GroupId::new("older"), 0.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.0, 0.0, 1.0),
            Vector3::try_new(0.0, 0.0, -1.0).unwrap(),
        ))
        .unwrap();
    let separation = MinimumFieldSeparation::try_new(1.0).unwrap();
    problem
        .add(YoungerThan::with_quadratic_penalty(
            SourceId::new("quadratic-age"),
            GroupId::new("younger"),
            GroupId::new("older"),
            separation,
            QuadraticPenalty::try_new(2.0).unwrap(),
        ))
        .unwrap();
    problem
        .add(YoungerThan::with_linear_violation_penalty(
            SourceId::new("linear-age"),
            GroupId::new("younger"),
            GroupId::new("older"),
            separation,
            LinearViolationPenalty::try_new(3.0).unwrap(),
        ))
        .unwrap();

    let success = problem.build().unwrap().fit().unwrap();
    let assessments = success.report().shared_level_set_relations();
    assert_eq!(assessments.len(), 2);
    assert_eq!(assessments[0].source_id().as_str(), "linear-age");
    assert_eq!(assessments[1].source_id().as_str(), "quadratic-age");
    for assessment in assessments {
        assert!(assessment.recovered_field_separation().abs() <= 1.0e-7);
        assert!((assessment.violation() - 1.0).abs() <= 1.0e-7);
        assert_eq!(assessment.active_state(), BoundActiveState::Active);
        assert!(assessment.loss().is_some());
    }
    assert_eq!(
        assessments[0].linear_violation_penalty().unwrap().weight(),
        3.0
    );
    assert!((assessments[0].loss().unwrap() - 3.0).abs() <= 1.0e-6);
    assert_eq!(assessments[1].quadratic_penalty().unwrap().weight(), 2.0);
    assert!((assessments[1].loss().unwrap() - 1.0).abs() <= 1.0e-6);
}

#[test]
fn duplicate_hard_age_relations_share_canonical_math_and_keep_both_sources() {
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    problem.add(horizon("older", "older/member", 0.0)).unwrap();
    problem
        .add(horizon("younger", "younger/member", 2.0))
        .unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("gauge"), GroupId::new("older"), 0.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.0, 0.0, 1.0),
            Vector3::try_new(0.0, 0.0, 1.0).unwrap(),
        ))
        .unwrap();
    let separation = MinimumFieldSeparation::try_new(1.0).unwrap();
    for source in ["age-b", "age-a"] {
        problem
            .add(YoungerThan::hard(
                SourceId::new(source),
                GroupId::new("younger"),
                GroupId::new("older"),
                separation,
            ))
            .unwrap();
    }

    let success = problem.build().unwrap().fit().unwrap();
    assert_eq!(success.report().shared_level_set_relations().len(), 2);
    assert_eq!(
        success
            .report()
            .problem_size()
            .affine_inequality_constraints(),
        1
    );
}

#[test]
fn accumulated_age_separation_conflicting_with_absolute_gauges_is_preflight_evidence() {
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    for (group, z) in [("a", 0.0), ("b", 1.0), ("c", 2.0)] {
        problem
            .add(horizon(group, &format!("{group}/member"), z))
            .unwrap();
    }
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("a-gauge"), GroupId::new("a"), 0.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("c-gauge"), GroupId::new("c"), 3.0)
                .unwrap(),
        )
        .unwrap();
    let separation = MinimumFieldSeparation::try_new(2.0).unwrap();
    problem
        .add(YoungerThan::hard(
            SourceId::new("a-to-b"),
            GroupId::new("b"),
            GroupId::new("a"),
            separation,
        ))
        .unwrap();
    problem
        .add(YoungerThan::hard(
            SourceId::new("b-to-c"),
            GroupId::new("c"),
            GroupId::new("b"),
            separation,
        ))
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    assert_eq!(
        failure.diagnosis(),
        georbf::diagnostics::ProblemDiagnosis::DirectInputConflict
    );
    let evidence = failure
        .report()
        .shared_level_set_relation_conflict()
        .unwrap();
    assert_eq!(
        evidence
            .source_ids()
            .iter()
            .map(SourceId::as_str)
            .collect::<Vec<_>>(),
        ["a-gauge", "a-to-b", "b-to-c", "c-gauge"]
    );
    assert_eq!(
        evidence
            .group_ids()
            .iter()
            .map(GroupId::as_str)
            .collect::<Vec<_>>(),
        ["a", "b", "c"]
    );
    let provenance = evidence.source_provenance();
    assert_eq!(provenance.len(), 4);
    assert_eq!(provenance[0].source_id().as_str(), "a-gauge");
    assert_eq!(
        provenance[0].semantic_role().as_str(),
        "additive-field-gauge/level-set"
    );
    assert_eq!(
        provenance[0]
            .group_ids()
            .iter()
            .map(GroupId::as_str)
            .collect::<Vec<_>>(),
        ["a"]
    );
    assert_eq!(provenance[1].source_id().as_str(), "a-to-b");
    assert_eq!(
        provenance[1]
            .group_ids()
            .iter()
            .map(GroupId::as_str)
            .collect::<Vec<_>>(),
        ["a", "b"]
    );
    assert_eq!(
        provenance[1].semantic_role().as_str(),
        "younger-than/minimum-field-separation"
    );
    assert_eq!(provenance[3].source_id().as_str(), "c-gauge");
    assert_eq!(
        provenance[3].semantic_role().as_str(),
        "additive-field-gauge/level-set"
    );
    assert!(failure.report().attempts().is_empty());
    let witness = failure.report().conflict_witness().unwrap();
    assert_eq!(witness.source_ids(), evidence.source_ids());
    assert_eq!(witness.canonical_residual(), 0.0);
    assert_eq!(witness.separation_margin(), 1.0);
    assert!(witness.provenance_verified());
    assert!(!witness.backend_invoked());
}

#[test]
fn absolute_field_value_conflict_retains_observation_and_member_roles() {
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    problem.add(horizon("a", "a-member", 0.0)).unwrap();
    problem.add(horizon("b", "b-member", 1.0)).unwrap();
    problem
        .add(
            FieldValueObservation::try_new(SourceId::new("a-absolute"), point(0.0, 0.0, 0.0), 0.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("b-gauge"), GroupId::new("b"), 1.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(YoungerThan::hard(
            SourceId::new("relation"),
            GroupId::new("b"),
            GroupId::new("a"),
            MinimumFieldSeparation::try_new(2.0).unwrap(),
        ))
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    assert_eq!(
        failure.diagnosis(),
        georbf::diagnostics::ProblemDiagnosis::DirectInputConflict
    );
    let evidence = failure
        .report()
        .shared_level_set_relation_conflict()
        .unwrap();
    let source = |source_id: &str| {
        evidence
            .source_provenance()
            .iter()
            .find(|source| source.source_id().as_str() == source_id)
            .unwrap()
    };
    assert_eq!(
        source("a-absolute").semantic_role().as_str(),
        "field-value-observation/value"
    );
    assert!(source("a-absolute").group_ids().is_empty());
    assert_eq!(
        source("a-member").semantic_role().as_str(),
        "shared-level-set/member/value"
    );
    assert_eq!(source("a-member").group_ids(), &[GroupId::new("a")]);
    assert_eq!(
        source("relation").semantic_role().as_str(),
        "younger-than/minimum-field-separation"
    );
    assert_eq!(
        source("relation").group_ids(),
        &[GroupId::new("a"), GroupId::new("b")]
    );
    assert!(!evidence.backend_invoked());
    assert!(failure.report().attempts().is_empty());
    let witness = failure.report().conflict_witness().unwrap();
    assert_eq!(witness.source_ids(), evidence.source_ids());
    assert_eq!(witness.canonical_residual(), 0.0);
    assert_eq!(witness.separation_margin(), 1.0);
}

#[test]
fn explicit_direction_alone_controls_age_semantics_across_ids_input_order_and_frame_metadata() {
    fn fit_case(
        older_group: &str,
        younger_group: &str,
        frame: InputCoordinateFrame,
        relation_first: bool,
    ) -> georbf::fit::FitSuccess {
        let mut problem = ProblemBuilder::new(frame, FieldUnitLabel::new("stratigraphic-unit"));
        problem
            .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
            .unwrap();
        let relation = || {
            OlderThan::hard(
                SourceId::new("age"),
                GroupId::new(older_group),
                GroupId::new(younger_group),
                MinimumFieldSeparation::try_new(1.0).unwrap(),
            )
        };
        if relation_first {
            problem.add(relation()).unwrap();
        }
        problem
            .add(horizon(younger_group, "younger/member", 2.0))
            .unwrap();
        problem
            .add(
                AdditiveFieldGauge::at_level_set(
                    SourceId::new("gauge"),
                    GroupId::new(older_group),
                    0.0,
                )
                .unwrap(),
            )
            .unwrap();
        problem
            .add(GradientObservation::new(
                SourceId::new("gradient"),
                point(0.0, 0.0, 1.0),
                Vector3::try_new(0.0, 0.0, 10.0).unwrap(),
            ))
            .unwrap();
        problem
            .add(horizon(older_group, "older/member", 0.0))
            .unwrap();
        problem
            .add(
                FieldValueObservation::try_new(
                    SourceId::new("absolute-younger"),
                    point(0.0, 0.0, 2.0),
                    20.0,
                )
                .unwrap(),
            )
            .unwrap();
        if !relation_first {
            problem.add(relation()).unwrap();
        }
        problem.build().unwrap().fit().unwrap()
    }

    let right = fit_case(
        "z-older",
        "a-younger",
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        true,
    );
    let left = fit_case(
        "a-older",
        "z-younger",
        InputCoordinateFrame::try_new(
            ["down", "west", "south"],
            Handedness::Left,
            LengthUnitLabel::new("m"),
        )
        .unwrap(),
        false,
    );
    for success in [&right, &left] {
        let query = success.model().evaluate(point(0.0, 0.0, 1.0)).unwrap();
        assert!((query.value() - 10.0).abs() <= 1.0e-6);
        assert!((query.gradient().components()[2] - 10.0).abs() <= 1.0e-6);
        let relation = &success.report().shared_level_set_relations()[0];
        assert_eq!(relation.kind(), SharedLevelSetRelationKind::OlderThan);
        assert!((relation.recovered_field_separation() - 20.0).abs() <= 1.0e-6);
        assert_eq!(
            relation.field_direction(),
            Some(StratigraphicFieldDirection::TowardYounger)
        );
    }
}

#[test]
fn general_hard_level_relation_infeasibility_requires_a_validated_backend_certificate() {
    let mut problem = builder();
    problem
        .set_stratigraphic_field_direction(StratigraphicFieldDirection::TowardYounger)
        .unwrap();
    problem.add(horizon("a", "a/member", 0.0)).unwrap();
    problem.add(horizon("b", "b/member", 1.0)).unwrap();
    problem
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.5, 0.0, 0.0),
            Vector3::try_new(0.0, 0.0, 0.0).unwrap(),
        ))
        .unwrap();
    problem
        .add(
            FieldValueBound::try_lower(SourceId::new("a-lower"), point(0.0, 0.0, 0.0), 0.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(
            FieldValueBound::try_upper(SourceId::new("b-upper"), point(0.0, 0.0, 1.0), 1.0)
                .unwrap(),
        )
        .unwrap();
    problem
        .add(YoungerThan::hard(
            SourceId::new("a-to-b"),
            GroupId::new("b"),
            GroupId::new("a"),
            MinimumFieldSeparation::try_new(2.0).unwrap(),
        ))
        .unwrap();

    let failure = problem.build().unwrap().fit().unwrap_err();
    assert_eq!(
        failure.diagnosis(),
        georbf::diagnostics::ProblemDiagnosis::InfeasibleProblem
    );
    let certificate = failure.report().infeasibility_certificate().unwrap();
    assert!(certificate.finite());
    assert!(certificate.backend_invoked());
    assert!(certificate.separation_margin() >= certificate.separation_limit());
    assert!(!failure.report().attempts().is_empty());
}
