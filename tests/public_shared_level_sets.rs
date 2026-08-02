use georbf::diagnostics::ProblemDiagnosis;
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::observation::{FieldValueObservation, GradientObservation};
use georbf::problem::{AddError, BuildError};
use georbf::relation::{
    AdditiveFieldGauge, GaugeError, GroupBuildError, GroupMemberAddError, HorizonBuilder,
    SharedLevelSetBuilder,
};
use georbf::{GroupId, Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the manufactured point is finite")
}

fn problem_builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .expect("the manufactured frame is valid"),
        FieldUnitLabel::new("stratigraphic-unit"),
    )
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 1.0e-9 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn complete_horizon_groups_are_added_atomically_after_forward_references() {
    let empty = SharedLevelSetBuilder::new(GroupId::new("empty"))
        .build()
        .expect_err("an empty shared level set cannot be constructed");
    assert_eq!(empty, GroupBuildError::EmptyGroup);

    let horizon_id = GroupId::new("horizon-a");
    let mut horizon = HorizonBuilder::new(horizon_id.clone());
    horizon
        .add_member(SourceId::new("horizon-a/0"), point(-1.0, 0.0, 0.5))
        .expect("the first member SourceId is unique");
    horizon
        .add_member(SourceId::new("horizon-a/1"), point(1.0, 0.0, -0.5))
        .expect("the second member SourceId is unique");

    let mut builder = problem_builder();
    builder
        .add(
            AdditiveFieldGauge::at_level_set(SourceId::new("gauge"), horizon_id.clone(), 4.0)
                .expect("the gauge value is finite"),
        )
        .expect("forward GroupId references are accepted by add");
    builder
        .add(horizon.build().expect("the horizon is non-empty"))
        .expect("the complete horizon is inserted in one atomic add");

    let snapshot = builder
        .build()
        .expect("the forward reference resolves by snapshot construction");
    assert_eq!(snapshot.horizon_count(), 1);
    assert_eq!(snapshot.shared_level_set_count(), 0);
    assert_eq!(snapshot.source_count(), 3);
}

#[test]
fn shift_invariant_horizons_report_the_unidentified_additive_gauge_before_solving() {
    let horizon_id = GroupId::new("shift-invariant-horizon");
    let mut horizon = HorizonBuilder::new(horizon_id.clone());
    horizon
        .add_member(SourceId::new("member-b"), point(1.0, 0.0, 0.0))
        .unwrap();
    horizon
        .add_member(SourceId::new("member-a"), point(-1.0, 0.0, 0.0))
        .unwrap();
    let mut builder = problem_builder();
    builder.add(horizon.build().unwrap()).unwrap();

    let failure = builder
        .build()
        .unwrap()
        .fit()
        .expect_err("relative shared-level relations do not choose an absolute representative");

    assert_eq!(
        failure.diagnosis(),
        ProblemDiagnosis::UnidentifiedAdditiveGauge
    );
    let evidence = failure
        .report()
        .unidentified_additive_gauge()
        .expect("the diagnosis retains stable semantic evidence");
    assert_eq!(evidence.group_ids(), &[horizon_id]);
    assert_eq!(
        evidence
            .source_ids()
            .iter()
            .map(SourceId::as_str)
            .collect::<Vec<_>>(),
        ["member-a", "member-b"]
    );
    assert!(!evidence.backend_invoked());
    assert!(failure.report().attempts().is_empty());
    assert_eq!(
        failure.report().problem_size().canonical_hard_equalities(),
        None
    );
    assert_eq!(failure.report().problem_size().kkt_dimension(), None);
}

#[test]
fn an_unreferenced_single_member_group_is_diagnosed_as_uninformative() {
    let group_id = GroupId::new("lonely-level-set");
    let mut group = SharedLevelSetBuilder::new(group_id.clone());
    group
        .add_member(SourceId::new("lonely-member"), point(0.0, 0.0, 0.0))
        .unwrap();
    let mut builder = problem_builder();
    builder.add(group.build().unwrap()).unwrap();

    let failure = builder
        .build()
        .unwrap()
        .fit()
        .expect_err("one disconnected member does not constrain its shared latent");

    assert_eq!(
        failure.diagnosis(),
        ProblemDiagnosis::UninformativeSharedLevelSet
    );
    let evidence = failure
        .report()
        .uninformative_shared_level_set()
        .expect("the diagnosis identifies the disconnected latent");
    assert_eq!(evidence.group_id(), &group_id);
    assert_eq!(evidence.member_source_id().as_str(), "lonely-member");
    assert!(!evidence.backend_invoked());
    assert!(failure.report().attempts().is_empty());
}

#[test]
fn repeated_locations_do_not_turn_a_multi_member_group_into_a_singleton() {
    let mut group = SharedLevelSetBuilder::new(GroupId::new("repeated-location-group"));
    group
        .add_member(SourceId::new("repeated/member-a"), point(0.0, 0.0, 0.0))
        .unwrap();
    group
        .add_member(SourceId::new("repeated/member-b"), point(0.0, 0.0, 0.0))
        .unwrap();
    let mut builder = problem_builder();
    builder.add(group.build().unwrap()).unwrap();

    let failure = builder.build().unwrap().fit().unwrap_err();
    assert_eq!(
        failure.diagnosis(),
        ProblemDiagnosis::UnidentifiedAdditiveGauge
    );
    assert!(failure.report().uninformative_shared_level_set().is_none());
}

#[test]
fn user_can_fit_a_manufactured_planar_horizon_and_recover_its_semantic_latent() {
    let horizon_id = GroupId::new("planar-horizon");
    let mut horizon = HorizonBuilder::new(horizon_id.clone());
    for (source, location) in [
        ("horizon/left", point(-2.0, 0.0, 1.0)),
        ("horizon/right", point(2.0, 0.0, -1.0)),
        ("horizon/north", point(0.0, 4.0, 1.0)),
    ] {
        horizon
            .add_member(SourceId::new(source), location)
            .expect("member SourceIds are unique");
    }
    let mut builder = problem_builder();
    builder.add(horizon.build().unwrap()).unwrap();
    builder
        .add(
            AdditiveFieldGauge::at_level_set(
                SourceId::new("gauge/planar-horizon"),
                horizon_id.clone(),
                3.0,
            )
            .unwrap(),
        )
        .unwrap();
    builder
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.5, -0.5, 0.25),
            Vector3::try_new(0.5, -0.25, 1.0).unwrap(),
        ))
        .unwrap();

    let success = builder
        .build()
        .unwrap()
        .fit()
        .expect("the gauge and complete gradient identify the manufactured plane");
    let model = success.model();
    assert_close(
        model
            .shared_level_value(&horizon_id)
            .expect("the semantic latent is addressable by GroupId"),
        3.0,
    );
    let query = point(0.75, -1.0, 2.0);
    let sample = model.evaluate(query).unwrap();
    assert_close(sample.value(), 3.0 + 0.5 * 0.75 - 0.25 * -1.0 + 2.0);
    for (actual, expected) in sample
        .gradient()
        .components()
        .into_iter()
        .zip([0.5, -0.25, 1.0])
    {
        assert_close(actual, expected);
    }

    let report = success.report();
    assert_eq!(report.problem_size().semantic_latents(), 1);
    assert_close(report.field_energy().unwrap(), 0.0);
    let recovered = report.shared_level_values();
    assert_eq!(recovered.len(), 1);
    assert_eq!(recovered[0].group_id(), &horizon_id);
    assert_close(recovered[0].value(), 3.0);
    assert_eq!(recovered[0].field_unit().as_str(), "stratigraphic-unit");
    assert_eq!(
        recovered[0]
            .member_source_ids()
            .iter()
            .map(SourceId::as_str)
            .collect::<Vec<_>>(),
        ["horizon/left", "horizon/north", "horizon/right"]
    );
    assert_eq!(
        report
            .hard_relations()
            .iter()
            .filter(|relation| relation.group_id() == Some(&horizon_id))
            .count(),
        4
    );
    let gauge = report
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id().as_str() == "gauge/planar-horizon")
        .expect("the convention retains independent source provenance");
    assert_eq!(gauge.group_id(), Some(&horizon_id));
    assert_eq!(
        gauge.semantic_role().as_str(),
        "additive-field-gauge/level-set"
    );
    assert!(
        report
            .hard_relations()
            .iter()
            .all(|relation| relation.residual().abs() <= relation.tolerance())
    );
}

fn fit_point_gauged_horizon(
    gauge_value: f64,
    reverse_input_order: bool,
) -> georbf::fit::FitSuccess {
    let horizon_id = GroupId::new("gauge-invariant-horizon");
    let mut horizon = HorizonBuilder::new(horizon_id);
    let mut members = [
        ("member-a", point(-2.0, 0.0, 1.0)),
        ("member-b", point(2.0, 0.0, -1.0)),
        ("member-c", point(0.0, 4.0, 1.0)),
    ];
    if reverse_input_order {
        members.reverse();
    }
    for (source, location) in members {
        horizon.add_member(SourceId::new(source), location).unwrap();
    }
    let horizon = horizon.build().unwrap();
    let gauge = AdditiveFieldGauge::at_point(
        SourceId::new("point-gauge"),
        point(0.0, 0.0, 0.0),
        gauge_value,
    )
    .unwrap();
    let gradient = GradientObservation::new(
        SourceId::new("gradient"),
        point(0.5, -0.5, 0.25),
        Vector3::try_new(0.5, -0.25, 1.0).unwrap(),
    );
    let mut builder = problem_builder();
    if reverse_input_order {
        builder.add(gauge).unwrap();
        builder.add(gradient).unwrap();
        builder.add(horizon).unwrap();
    } else {
        builder.add(horizon).unwrap();
        builder.add(gradient).unwrap();
        builder.add(gauge).unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

fn fit_horizon_with_point_gauges(
    first_value: f64,
    second_value: Option<f64>,
) -> Result<georbf::fit::FitSuccess, georbf::fit::FitFailure> {
    let mut horizon = HorizonBuilder::new(GroupId::new("multi-gauge-horizon"));
    for (source, location) in [
        ("multi/member-a", point(-2.0, 0.0, 1.0)),
        ("multi/member-b", point(2.0, 0.0, -1.0)),
        ("multi/member-c", point(0.0, 4.0, 1.0)),
    ] {
        horizon.add_member(SourceId::new(source), location).unwrap();
    }
    let mut builder = problem_builder();
    builder.add(horizon.build().unwrap()).unwrap();
    builder
        .add(GradientObservation::new(
            SourceId::new("multi/gradient"),
            point(0.5, -0.5, 0.25),
            Vector3::try_new(0.5, -0.25, 1.0).unwrap(),
        ))
        .unwrap();
    builder
        .add(
            AdditiveFieldGauge::at_point(
                SourceId::new("multi/gauge-a"),
                point(0.0, 0.0, 0.0),
                first_value,
            )
            .unwrap(),
        )
        .unwrap();
    if let Some(value) = second_value {
        builder
            .add(
                AdditiveFieldGauge::at_point(
                    SourceId::new("multi/gauge-b"),
                    point(2.0, 0.0, 0.0),
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder.build().unwrap().fit()
}

#[test]
fn additive_gauge_and_input_reordering_preserve_canonical_observables() {
    let baseline = fit_point_gauged_horizon(3.0, false);
    let reordered = fit_point_gauged_horizon(3.0, true);
    let shifted = fit_point_gauged_horizon(10.0, true);
    let group_id = GroupId::new("gauge-invariant-horizon");
    let query = point(0.75, -1.0, 2.0);

    let baseline_sample = baseline.model().evaluate(query).unwrap();
    let reordered_sample = reordered.model().evaluate(query).unwrap();
    let shifted_sample = shifted.model().evaluate(query).unwrap();
    assert_close(reordered_sample.value(), baseline_sample.value());
    assert_close(shifted_sample.value() - baseline_sample.value(), 7.0);
    for ((baseline_component, reordered_component), shifted_component) in baseline_sample
        .gradient()
        .components()
        .into_iter()
        .zip(reordered_sample.gradient().components())
        .zip(shifted_sample.gradient().components())
    {
        assert_close(reordered_component, baseline_component);
        assert_close(shifted_component, baseline_component);
    }
    assert_close(
        reordered.model().shared_level_value(&group_id).unwrap(),
        baseline.model().shared_level_value(&group_id).unwrap(),
    );
    assert_close(
        shifted.model().shared_level_value(&group_id).unwrap()
            - baseline.model().shared_level_value(&group_id).unwrap(),
        7.0,
    );
    assert_close(
        reordered.report().field_energy().unwrap(),
        baseline.report().field_energy().unwrap(),
    );
    assert_close(
        shifted.report().field_energy().unwrap(),
        baseline.report().field_energy().unwrap(),
    );
    assert_eq!(
        reordered.report().shared_level_values(),
        baseline.report().shared_level_values()
    );
}

#[test]
fn additional_gauges_verify_the_additive_representative_without_changing_geometry() {
    let baseline = fit_horizon_with_point_gauges(3.0, None).unwrap();
    let compatible = fit_horizon_with_point_gauges(3.0, Some(4.0)).unwrap();
    let shifted = fit_horizon_with_point_gauges(10.0, Some(11.0)).unwrap();
    let query = point(0.75, -1.0, 2.0);
    let baseline_sample = baseline.model().evaluate(query).unwrap();
    let compatible_sample = compatible.model().evaluate(query).unwrap();
    let shifted_sample = shifted.model().evaluate(query).unwrap();

    for ((baseline_component, compatible_component), shifted_component) in baseline_sample
        .gradient()
        .components()
        .into_iter()
        .zip(compatible_sample.gradient().components())
        .zip(shifted_sample.gradient().components())
    {
        assert_close(compatible_component, baseline_component);
        assert_close(shifted_component, baseline_component);
    }
    assert_close(
        compatible.report().field_energy().unwrap(),
        baseline.report().field_energy().unwrap(),
    );
    assert_close(
        shifted.report().field_energy().unwrap(),
        baseline.report().field_energy().unwrap(),
    );
    assert_close(compatible_sample.value(), baseline_sample.value());
    assert_close(shifted_sample.value() - baseline_sample.value(), 7.0);
    assert_eq!(
        compatible
            .report()
            .hard_relations()
            .iter()
            .filter(|relation| relation.source_id().as_str().starts_with("multi/gauge-"))
            .count(),
        2
    );
    let primary_gauge = compatible
        .report()
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id().as_str() == "multi/gauge-a")
        .unwrap();
    let verification_gauge = compatible
        .report()
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id().as_str() == "multi/gauge-b")
        .unwrap();
    assert!(primary_gauge.scaled_kkt_tolerance().is_some());
    assert!(verification_gauge.scaled_kkt_tolerance().is_none());

    let failure = fit_horizon_with_point_gauges(3.0, Some(4.5))
        .expect_err("an incompatible secondary convention must reject during recovery");
    assert_eq!(
        failure.diagnosis(),
        ProblemDiagnosis::RecoveryVerificationFailure
    );
    assert!(failure.report().recovery_verification().is_some());
    let conflicting = failure
        .report()
        .hard_relations()
        .iter()
        .find(|relation| relation.source_id().as_str() == "multi/gauge-b")
        .expect("the failed verification retains secondary-gauge provenance");
    assert!(conflicting.residual().abs() > conflicting.tolerance());
}

#[test]
fn absolute_observations_close_group_cycles_without_redundant_kkt_rows() {
    let group_id = GroupId::new("absolutely-observed-horizon");
    let left = point(-1.0, 0.0, 0.5);
    let right = point(1.0, 0.0, -0.5);
    let mut horizon = HorizonBuilder::new(group_id.clone());
    horizon
        .add_member(SourceId::new("absolute/member-left"), left)
        .unwrap();
    horizon
        .add_member(SourceId::new("absolute/member-right"), right)
        .unwrap();
    let mut builder = problem_builder();
    builder.add(horizon.build().unwrap()).unwrap();
    builder
        .add(
            FieldValueObservation::try_new(SourceId::new("absolute/value-left"), left, 3.0)
                .unwrap(),
        )
        .unwrap();
    builder
        .add(
            FieldValueObservation::try_new(SourceId::new("absolute/value-right"), right, 3.0)
                .unwrap(),
        )
        .unwrap();
    builder
        .add(GradientObservation::new(
            SourceId::new("absolute/gradient"),
            point(0.0, 1.0, 0.0),
            Vector3::try_new(0.5, -0.25, 1.0).unwrap(),
        ))
        .unwrap();

    let success = builder
        .build()
        .unwrap()
        .fit()
        .expect("the compatible absolute/member cycle should verify without rank loss");
    assert_close(success.model().shared_level_value(&group_id).unwrap(), 3.0);
    assert!(
        success
            .report()
            .hard_relations()
            .iter()
            .all(|relation| relation.residual().abs() <= relation.tolerance())
    );
    assert_eq!(success.report().problem_size().input_observations(), 3);
    assert_eq!(
        success.report().problem_size().canonical_hard_equalities(),
        Some(7)
    );
    assert_eq!(
        success.report().problem_size().equality_constraints(),
        Some(10)
    );

    let mut conflicting_horizon = HorizonBuilder::new(group_id);
    conflicting_horizon
        .add_member(SourceId::new("conflict/member-left"), left)
        .unwrap();
    conflicting_horizon
        .add_member(SourceId::new("conflict/member-right"), right)
        .unwrap();
    let mut conflicting = problem_builder();
    conflicting
        .add(conflicting_horizon.build().unwrap())
        .unwrap();
    conflicting
        .add(
            FieldValueObservation::try_new(SourceId::new("conflict/value-left"), left, 3.0)
                .unwrap(),
        )
        .unwrap();
    conflicting
        .add(
            FieldValueObservation::try_new(SourceId::new("conflict/value-right"), right, 4.0)
                .unwrap(),
        )
        .unwrap();
    let failure = conflicting
        .build()
        .unwrap()
        .fit()
        .expect_err("a graph-provable shared-level contradiction must fail pre-backend");
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().direct_input_conflict().is_some());
    assert!(failure.report().attempts().is_empty());
}

#[test]
fn a_level_set_gauge_makes_a_single_member_group_informative() {
    let group_id = GroupId::new("referenced-singleton");
    let mut group = SharedLevelSetBuilder::new(group_id.clone());
    group
        .add_member(SourceId::new("singleton"), point(0.0, 0.0, 0.0))
        .unwrap();
    let mut builder = problem_builder();
    builder.add(group.build().unwrap()).unwrap();
    builder
        .add(
            AdditiveFieldGauge::at_level_set(
                SourceId::new("singleton-gauge"),
                group_id.clone(),
                3.0,
            )
            .unwrap(),
        )
        .unwrap();
    for (index, location) in [point(1.0, 0.0, 0.0), point(0.0, 1.0, 1.0)]
        .into_iter()
        .enumerate()
    {
        builder
            .add(GradientObservation::new(
                SourceId::new(format!("singleton-gradient-{index}")),
                location,
                Vector3::try_new(0.5, -0.25, 1.0).unwrap(),
            ))
            .unwrap();
    }

    let success = builder
        .build()
        .unwrap()
        .fit()
        .expect("the explicit group reference makes the singleton latent informative");
    assert_close(success.model().shared_level_value(&group_id).unwrap(), 3.0);
}

#[test]
fn group_and_gauge_validation_is_atomic_and_repairable() {
    assert_eq!(
        AdditiveFieldGauge::at_point(SourceId::new("bad-gauge"), point(0.0, 0.0, 0.0), f64::NAN,)
            .unwrap_err(),
        GaugeError::NonFiniteFieldValue
    );

    let mut group = SharedLevelSetBuilder::new(GroupId::new("target"));
    group
        .add_member(SourceId::new("member"), point(0.0, 0.0, 0.0))
        .unwrap();
    assert!(matches!(
        group.add_member(SourceId::new("member"), point(1.0, 0.0, 0.0)),
        Err(GroupMemberAddError::DuplicateSourceId { .. })
    ));
    let group = group
        .build()
        .expect("the rejected member did not mutate the draft");

    let mut builder = problem_builder();
    builder
        .add(
            AdditiveFieldGauge::at_level_set(
                SourceId::new("forward-gauge"),
                GroupId::new("target"),
                0.0,
            )
            .unwrap(),
        )
        .unwrap();
    let failure = builder
        .build()
        .expect_err("a still-dangling forward reference is a build error");
    assert!(matches!(
        failure.errors(),
        [BuildError::UnknownGroupReference { group_id, .. }] if group_id.as_str() == "target"
    ));
    let mut repaired = failure.into_builder();
    repaired.add(group.clone()).unwrap();

    let mut duplicate = SharedLevelSetBuilder::new(GroupId::new("target"));
    duplicate
        .add_member(SourceId::new("other-member"), point(1.0, 0.0, 0.0))
        .unwrap();
    assert!(matches!(
        repaired.add(duplicate.build().unwrap()),
        Err(AddError::DuplicateGroupId { .. })
    ));
    assert_eq!(repaired.build().unwrap().source_count(), 2);
}

#[test]
fn overlapping_member_cycles_do_not_merge_semantic_latent_identity() {
    let mut builder = problem_builder();
    for group_name in ["group-a", "group-b"] {
        let mut group = SharedLevelSetBuilder::new(GroupId::new(group_name));
        let shared_zero = if group_name == "group-b" { -0.0 } else { 0.0 };
        group
            .add_member(
                SourceId::new(format!("{group_name}/left")),
                point(-1.0, shared_zero, 0.5),
            )
            .unwrap();
        group
            .add_member(
                SourceId::new(format!("{group_name}/right")),
                point(1.0, shared_zero, -0.5),
            )
            .unwrap();
        if group_name == "group-b" {
            group
                .add_member(SourceId::new("group-b/north"), point(0.0, 2.0, 0.5))
                .unwrap();
        }
        builder.add(group.build().unwrap()).unwrap();
    }
    builder
        .add(
            AdditiveFieldGauge::at_point(SourceId::new("shared-gauge"), point(0.0, 0.0, 0.0), 3.0)
                .unwrap(),
        )
        .unwrap();
    builder
        .add(GradientObservation::new(
            SourceId::new("shared-gradient"),
            point(0.0, 1.0, 0.0),
            Vector3::try_new(0.5, -0.25, 1.0).unwrap(),
        ))
        .unwrap();

    let success = builder.build().unwrap().fit().unwrap();
    assert_eq!(success.report().problem_size().semantic_latents(), 2);
    assert_eq!(
        success.report().problem_size().canonical_hard_equalities(),
        Some(9)
    );
    assert_close(
        success
            .model()
            .shared_level_value(&GroupId::new("group-a"))
            .unwrap(),
        3.0,
    );
    assert_close(
        success
            .model()
            .shared_level_value(&GroupId::new("group-b"))
            .unwrap(),
        3.0,
    );
}
