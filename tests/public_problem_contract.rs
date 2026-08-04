use std::num::NonZeroUsize;
use std::thread;

use georbf::geometry::{
    FieldUnitLabel, GeometryError, GlobalAnisotropyMetric, Handedness, InputCoordinateFrame,
    LengthUnitLabel, Point3, Vector3,
};
use georbf::observation::{
    FieldValueObservation, GradientObservation, ObservationError, TangentDirectionObservation,
};
use georbf::problem::{AddError, BuildError, FitConfiguration, ThreadBudget};
use georbf::relation::{
    AdditiveFieldGauge, GaugeError, GroupBuildError, GroupMemberAddError, HorizonBuilder,
    SharedLevelSetBuilder,
};
use georbf::{GroupId, ProblemBuilder, SourceId};
use proptest::prelude::*;

fn frame() -> InputCoordinateFrame {
    InputCoordinateFrame::try_new(
        ["east", "north", "elevation"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )
    .expect("the fixture frame is valid")
}

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the fixture point is finite")
}

fn affine_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    2.0 + 0.5 * x - 1.25 * y + 0.75 * z
}

#[test]
fn property_profile_case_override_is_positive_and_fail_closed() {
    assert_eq!(configured_property_cases(128, None), Ok(128));
    assert_eq!(configured_property_cases(128, Some("2500")), Ok(2_500));
    assert!(configured_property_cases(128, Some("0")).is_err());
    assert!(configured_property_cases(128, Some("not-a-number")).is_err());
}

fn configured_property_cases(default: u32, override_value: Option<&str>) -> Result<u32, &str> {
    let Some(value) = override_value else {
        return Ok(default);
    };
    value
        .parse::<u32>()
        .ok()
        .filter(|cases| *cases > 0)
        .ok_or("PROPTEST_CASES must be a positive u32")
}

fn property_cases(default: u32) -> u32 {
    match std::env::var("PROPTEST_CASES") {
        Ok(value) => configured_property_cases(default, Some(&value))
            .expect("PROPTEST_CASES must be a positive u32"),
        Err(std::env::VarError::NotPresent) => default,
        Err(std::env::VarError::NotUnicode(_)) => {
            panic!("PROPTEST_CASES must contain valid Unicode")
        }
    }
}

fn add_affine_input(builder: &mut ProblemBuilder, input_index: usize) {
    let value_locations = [
        point(-1.0, -1.0, -1.0),
        point(1.0, -1.0, -1.0),
        point(-1.0, 1.0, -1.0),
        point(-1.0, -1.0, 1.0),
        point(1.0, 1.0, 0.5),
    ];
    let gradient_locations = [
        point(0.25, -0.5, 0.75),
        point(-0.75, 0.25, 0.5),
        point(0.5, 0.75, -0.25),
    ];
    match input_index {
        0..=4 => {
            let location = value_locations[input_index];
            builder
                .add(
                    FieldValueObservation::try_new(
                        SourceId::new(format!("value-{input_index}")),
                        location,
                        affine_value(location),
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        5 => {
            let location = value_locations[0];
            builder
                .add(
                    FieldValueObservation::try_new(
                        SourceId::new("value-alias"),
                        location,
                        affine_value(location),
                    )
                    .unwrap(),
                )
                .unwrap();
        }
        6..=8 => {
            let gradient_index = input_index - 6;
            builder
                .add(GradientObservation::new(
                    SourceId::new(format!("gradient-{gradient_index}")),
                    gradient_locations[gradient_index],
                    Vector3::try_new(0.5, -1.25, 0.75).unwrap(),
                ))
                .unwrap();
        }
        _ => panic!("the affine fixture has exactly nine inputs"),
    }
}

fn affine_snapshot(order: &[usize]) -> georbf::ProblemSnapshot {
    let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field"));
    for input_index in order {
        add_affine_input(&mut builder, *input_index);
    }
    builder.build().expect("the affine fixture is valid")
}

macro_rules! assert_group_draft_contract {
    ($builder:ty, $empty_id:literal, $group_id:literal, $member_count:expr, $duplicate_index:expr) => {{
        let empty = <$builder>::new(GroupId::new($empty_id));
        prop_assert_eq!(empty.build().unwrap_err(), GroupBuildError::EmptyGroup);

        let mut draft = <$builder>::new(GroupId::new($group_id));
        for index in 0..$member_count {
            draft
                .add_member(
                    SourceId::new(format!("member-{index:02}")),
                    point(index as f64, 0.0, 0.0),
                )
                .unwrap();
        }
        let duplicate_id = SourceId::new(format!("member-{:02}", $duplicate_index));
        prop_assert_eq!(
            draft.add_member(duplicate_id.clone(), point(99.0, 0.0, 0.0)),
            Err(GroupMemberAddError::DuplicateSourceId {
                source_id: duplicate_id,
            })
        );
        prop_assert_eq!(draft.build().unwrap().members().len(), $member_count);
    }};
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(property_cases(256)))]

    #[test]
    fn checked_leaf_constructors_accept_exactly_their_finite_domains(
        components in any::<[u64; 3]>(),
        value_bits in any::<u64>(),
    ) {
        let [x, y, z] = components.map(f64::from_bits);
        let value = f64::from_bits(value_bits);
        let all_components_are_finite = [x, y, z].into_iter().all(f64::is_finite);

        prop_assert_eq!(Point3::try_new(x, y, z).is_ok(), all_components_are_finite);
        prop_assert_eq!(Vector3::try_new(x, y, z).is_ok(), all_components_are_finite);

        let location = point(0.0, 0.0, 0.0);
        let field_value = FieldValueObservation::try_new(
            SourceId::new("field-value"),
            location,
            value,
        );
        let gauge = AdditiveFieldGauge::at_point(
            SourceId::new("gauge"),
            location,
            value,
        );
        prop_assert_eq!(field_value.is_ok(), value.is_finite());
        prop_assert_eq!(gauge.is_ok(), value.is_finite());
        if !value.is_finite() {
            prop_assert_eq!(field_value.unwrap_err(), ObservationError::NonFiniteFieldValue);
            prop_assert_eq!(gauge.unwrap_err(), GaugeError::NonFiniteFieldValue);
        }
    }

    #[test]
    fn tangent_construction_never_panics_for_finite_vector_bit_patterns(
        components in any::<[u64; 3]>(),
    ) {
        let [x, y, z] = components.map(f64::from_bits);
        let Ok(direction) = Vector3::try_new(x, y, z) else {
            return Ok(());
        };
        let result = TangentDirectionObservation::try_new(
            SourceId::new("tangent"),
            point(0.0, 0.0, 0.0),
            direction,
        );
        if [x, y, z].into_iter().all(|component| component == 0.0) {
            prop_assert_eq!(result.unwrap_err(), ObservationError::ZeroTangentDirection);
        } else {
            let unit = result.expect("every finite nonzero vector has an axial representative");
            let components = unit.direction().components();
            prop_assert!(components.into_iter().all(f64::is_finite));
            let norm = components
                .into_iter()
                .map(|component| component * component)
                .sum::<f64>()
                .sqrt();
            prop_assert!((norm - 1.0).abs() <= 8.0 * f64::EPSILON);
            let first_nonzero = components.into_iter().find(|component| *component != 0.0);
            prop_assert!(first_nonzero.is_none_or(f64::is_sign_positive));
        }
    }

    #[test]
    fn frame_validation_is_a_total_checked_operation(labels in prop::array::uniform3(0_u8..4)) {
        let labels = labels.map(|label| match label {
            0 => String::new(),
            value => format!("axis-{value}"),
        });
        let expected_valid = labels.iter().all(|label| !label.is_empty())
            && labels[0] != labels[1]
            && labels[0] != labels[2]
            && labels[1] != labels[2];
        let result = InputCoordinateFrame::try_new(
            labels,
            Handedness::Right,
            LengthUnitLabel::new("metre"),
        );
        prop_assert_eq!(result.is_ok(), expected_valid);
        if let Err(error) = result {
            let is_structured_frame_error = matches!(
                error,
                GeometryError::EmptyAxisLabel { .. }
                    | GeometryError::DuplicateAxisLabel { .. }
            );
            prop_assert!(is_structured_frame_error);
        }
    }

    #[test]
    fn anisotropy_validation_never_repairs_or_panics_for_arbitrary_matrices(
        matrix_bits in any::<[[u64; 3]; 3]>(),
    ) {
        let matrix = matrix_bits.map(|row| row.map(f64::from_bits));
        if let Ok(metric) = GlobalAnisotropyMetric::try_from_matrix(matrix) {
            prop_assert_eq!(metric.matrix(), matrix);
            prop_assert!(matrix.into_iter().flatten().all(f64::is_finite));
            prop_assert_eq!(matrix[0][1], matrix[1][0]);
            prop_assert_eq!(matrix[0][2], matrix[2][0]);
            prop_assert_eq!(matrix[1][2], matrix[2][1]);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(property_cases(128)))]

    #[test]
    fn duplicate_source_rejection_leaves_the_builder_unchanged(
        unique_count in 1_usize..16,
        duplicate_choice in any::<usize>(),
    ) {
        let duplicate_index = duplicate_choice % unique_count;
        let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field"));
        for index in 0..unique_count {
            builder
                .add(GradientObservation::new(
                    SourceId::new(format!("source-{index:02}")),
                    point(index as f64, 0.0, 0.0),
                    Vector3::try_new(1.0, 0.0, 0.0).unwrap(),
                ))
                .unwrap();
        }

        let duplicate_id = SourceId::new(format!("source-{duplicate_index:02}"));
        let rejected = builder.add(
            FieldValueObservation::try_new(duplicate_id.clone(), point(99.0, 0.0, 0.0), 4.0)
                .unwrap(),
        );
        prop_assert_eq!(
            rejected,
            Err(AddError::DuplicateSourceId {
                source_id: duplicate_id,
            })
        );

        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new("post-rejection"),
                    point(-1.0, 0.0, 0.0),
                    3.0,
                )
                .unwrap(),
            )
            .expect("a rejected add cannot reserve or append any state");
        let snapshot = builder.build().unwrap();
        prop_assert_eq!(snapshot.observation_count(), unique_count + 1);
        prop_assert_eq!(snapshot.source_count(), unique_count + 1);
    }

    #[test]
    fn dangling_references_aggregate_stably_and_the_builder_can_be_repaired(
        priorities in prop::collection::vec(any::<u16>(), 1..8),
        requested_threads in 2_usize..9,
    ) {
        let count = priorities.len();
        let mut order = (0..count).collect::<Vec<_>>();
        order.sort_by_key(|index| (priorities[*index], *index));

        let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field"));
        builder.set_fit_configuration(
            FitConfiguration::default().with_thread_budget(ThreadBudget::Exact(
                NonZeroUsize::new(requested_threads).unwrap(),
            )),
        );
        for index in order {
            builder
                .add(
                    AdditiveFieldGauge::at_level_set(
                        SourceId::new(format!("gauge-{index:02}")),
                        GroupId::new(format!("group-{index:02}")),
                        index as f64,
                    )
                    .unwrap(),
                )
                .unwrap();
        }

        let failure = builder
            .build()
            .expect_err("all missing groups and the unsupported resource request are reported");
        let mut expected = (0..count)
            .map(|index| BuildError::UnknownGroupReference {
                source_id: SourceId::new(format!("gauge-{index:02}")),
                group_id: GroupId::new(format!("group-{index:02}")),
            })
            .collect::<Vec<_>>();
        expected.push(BuildError::UnsupportedThreadBudget {
            requested: requested_threads,
        });
        prop_assert_eq!(failure.errors(), expected.as_slice());

        let mut repaired = failure.into_builder();
        repaired.set_fit_configuration(FitConfiguration::default());
        for index in (0..count).rev() {
            let mut group = SharedLevelSetBuilder::new(GroupId::new(format!("group-{index:02}")));
            group
                .add_member(
                    SourceId::new(format!("member-{index:02}")),
                    point(index as f64, 1.0, 0.0),
                )
                .unwrap();
            repaired.add(group.build().unwrap()).unwrap();
        }
        let snapshot = repaired
            .build()
            .expect("repair preserves forward references and all original gauges");
        prop_assert_eq!(snapshot.shared_level_set_count(), count);
        prop_assert_eq!(snapshot.source_count(), count * 2);
    }

    #[test]
    fn incomplete_groups_and_duplicate_members_are_checked_atomically(
        horizon in any::<bool>(),
        member_count in 1_usize..12,
        duplicate_choice in any::<usize>(),
    ) {
        let duplicate_index = duplicate_choice % member_count;
        if horizon {
            assert_group_draft_contract!(
                HorizonBuilder,
                "empty-horizon",
                "horizon",
                member_count,
                duplicate_index
            );
        } else {
            assert_group_draft_contract!(
                SharedLevelSetBuilder,
                "empty-level-set",
                "level-set",
                member_count,
                duplicate_index
            );
        }
    }

    #[test]
    fn completed_groups_reject_duplicate_group_and_source_ids_without_partial_mutation(
        member_count in 1_usize..10,
        collision_choice in any::<usize>(),
    ) {
        let collision_index = collision_choice % member_count;
        let mut accepted = SharedLevelSetBuilder::new(GroupId::new("accepted-group"));
        for index in 0..member_count {
            accepted
                .add_member(
                    SourceId::new(format!("accepted-member-{index:02}")),
                    point(index as f64, 0.0, 0.0),
                )
                .unwrap();
        }
        let mut builder = ProblemBuilder::new(frame(), FieldUnitLabel::new("field"));
        builder.add(accepted.build().unwrap()).unwrap();

        let mut duplicate_group = SharedLevelSetBuilder::new(GroupId::new("accepted-group"));
        duplicate_group
            .add_member(SourceId::new("new-member"), point(20.0, 0.0, 0.0))
            .unwrap();
        prop_assert_eq!(
            builder.add(duplicate_group.build().unwrap()),
            Err(AddError::DuplicateGroupId {
                group_id: GroupId::new("accepted-group"),
            })
        );

        let colliding_source = SourceId::new(format!("accepted-member-{collision_index:02}"));
        let mut duplicate_source = SharedLevelSetBuilder::new(GroupId::new("new-group"));
        duplicate_source
            .add_member(colliding_source.clone(), point(21.0, 0.0, 0.0))
            .unwrap();
        duplicate_source
            .add_member(SourceId::new("also-not-added"), point(22.0, 0.0, 0.0))
            .unwrap();
        prop_assert_eq!(
            builder.add(duplicate_source.build().unwrap()),
            Err(AddError::DuplicateSourceId {
                source_id: colliding_source,
            })
        );

        let snapshot = builder.build().unwrap();
        prop_assert_eq!(snapshot.shared_level_set_count(), 1);
        prop_assert_eq!(snapshot.source_count(), member_count);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(property_cases(32)))]

    #[test]
    fn input_permutation_preserves_build_fit_and_exact_duplicate_identity(
        priorities in prop::collection::vec(any::<u16>(), 9),
    ) {
        let baseline_order = (0..9).collect::<Vec<_>>();
        let mut permuted_order = baseline_order.clone();
        permuted_order.sort_by_key(|index| (priorities[*index], *index));

        let baseline_snapshot = affine_snapshot(&baseline_order);
        let permuted_snapshot = affine_snapshot(&permuted_order);
        prop_assert_eq!(
            permuted_snapshot.input_coordinate_frame(),
            baseline_snapshot.input_coordinate_frame()
        );
        prop_assert_eq!(permuted_snapshot.field_unit(), baseline_snapshot.field_unit());
        prop_assert_eq!(permuted_snapshot.observation_count(), 9);
        prop_assert_eq!(permuted_snapshot.source_count(), 9);

        let baseline = baseline_snapshot.fit().expect("the affine baseline fits");
        let permuted = permuted_snapshot.fit().expect("the permuted affine problem fits");
        prop_assert_eq!(
            permuted.report().problem_size(),
            baseline.report().problem_size()
        );
        prop_assert_eq!(
            permuted.report().hard_relations(),
            baseline.report().hard_relations()
        );
        prop_assert_eq!(
            permuted.report().shared_level_values(),
            baseline.report().shared_level_values()
        );
        prop_assert_eq!(
            permuted.report().field_energy(),
            baseline.report().field_energy()
        );
        prop_assert_eq!(
            permuted.report().total_objective(),
            baseline.report().total_objective()
        );
        let query = point(0.2, -0.3, 0.4);
        prop_assert_eq!(
            permuted.model().evaluate(query).unwrap(),
            baseline.model().evaluate(query).unwrap()
        );
        prop_assert!(
            permuted
                .report()
                .canonical_acceptance()
                .is_some_and(|acceptance| acceptance.provenance_verified())
        );

        let report = permuted.report();
        prop_assert_eq!(report.problem_size().scalar_hard_relations(), 15);
        prop_assert_eq!(report.problem_size().canonical_hard_equalities(), Some(14));
        prop_assert_eq!(
            report
                .hard_relations()
                .iter()
                .filter(|relation| matches!(
                    relation.source_id().as_str(),
                    "value-0" | "value-alias"
                ))
                .count(),
            2
        );
    }
}

#[test]
fn snapshots_are_send_sync_and_preserve_declared_problem_metadata() {
    fn assert_send_sync_clone<T: Send + Sync + Clone>() {}
    assert_send_sync_clone::<georbf::ProblemSnapshot>();

    let input_frame = frame();
    let field_unit = FieldUnitLabel::new("stratigraphic-coordinate");
    let mut builder = ProblemBuilder::new(input_frame.clone(), field_unit.clone());
    builder
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.0, 0.0, 0.0),
            Vector3::try_new(1.0, 2.0, 3.0).unwrap(),
        ))
        .unwrap();
    let snapshot = builder.build().unwrap();

    assert_eq!(snapshot.input_coordinate_frame(), &input_frame);
    assert_eq!(snapshot.field_unit(), &field_unit);
    assert_eq!(snapshot.observation_count(), 1);
    assert_eq!(snapshot.source_count(), 1);
}

#[test]
fn cloned_snapshots_support_concurrent_readonly_fits() {
    let snapshot = affine_snapshot(&(0..9).collect::<Vec<_>>());
    let expected = snapshot.fit().unwrap();
    let expected_relations = expected.report().hard_relations().to_vec();
    let query = point(0.2, -0.3, 0.4);
    let expected_sample = expected.model().evaluate(query).unwrap();

    let handles = (0..4)
        .map(|_| {
            let snapshot = snapshot.clone();
            thread::spawn(move || {
                let success = snapshot.fit().expect("a read-only cloned snapshot fits");
                (
                    success.report().hard_relations().to_vec(),
                    success.model().evaluate(query).unwrap(),
                    snapshot.input_coordinate_frame().clone(),
                    snapshot.field_unit().clone(),
                )
            })
        })
        .collect::<Vec<_>>();

    for handle in handles {
        let (relations, sample, input_frame, field_unit) =
            handle.join().expect("the snapshot read thread succeeds");
        assert_eq!(relations, expected_relations);
        assert_eq!(sample, expected_sample);
        assert_eq!(input_frame, frame());
        assert_eq!(field_unit, FieldUnitLabel::new("field"));
    }
}
