use georbf::diagnostics::{
    ProblemDiagnosis, RankDecision, RankDeficiencyConcept, RankEvidenceDomain,
    SolveAttemptTermination,
};
use georbf::geometry::{FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel};
use georbf::observation::{FieldValueObservation, GradientObservation};
use georbf::relation::SharedLevelSetBuilder;
use georbf::{GroupId, Point3, ProblemBuilder, SourceId, Vector3};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the diagnostic fixture point is finite")
}

fn problem_builder() -> ProblemBuilder {
    ProblemBuilder::new(
        InputCoordinateFrame::try_new(
            ["east", "north", "elevation"],
            Handedness::Right,
            LengthUnitLabel::new("m"),
        )
        .expect("the diagnostic fixture frame is valid"),
        FieldUnitLabel::new("stratigraphic-unit"),
    )
}

fn singleton(
    group_id: &str,
    source_id: &str,
    location: Point3,
) -> georbf::relation::SharedLevelSet {
    let mut group = SharedLevelSetBuilder::new(GroupId::new(group_id));
    group
        .add_member(SourceId::new(source_id), location)
        .expect("the singleton source is unique");
    group.build().expect("the singleton is nonempty")
}

fn mixed_preflight_failure(reverse_input_order: bool) -> georbf::fit::FitFailure {
    let group_a = singleton("a-group", "z-member", point(-1.0, 0.0, 0.0));
    let group_z = singleton("z-group", "a-member", point(1.0, 0.0, 0.0));
    let value_a =
        FieldValueObservation::try_new(SourceId::new("conflict-a"), point(0.0, 0.0, 0.0), 1.0)
            .unwrap();
    let value_b =
        FieldValueObservation::try_new(SourceId::new("conflict-b"), point(0.0, 0.0, 0.0), 2.0)
            .unwrap();

    let mut builder = problem_builder();
    if reverse_input_order {
        builder.add(value_b).unwrap();
        builder.add(group_z).unwrap();
        builder.add(value_a).unwrap();
        builder.add(group_a).unwrap();
    } else {
        builder.add(group_a).unwrap();
        builder.add(value_a).unwrap();
        builder.add(group_z).unwrap();
        builder.add(value_b).unwrap();
    }
    builder
        .build()
        .unwrap()
        .fit()
        .expect_err("unresolved semantics and an exact conflict must fail before solving")
}

#[test]
fn preflight_retains_all_evidence_and_selects_primary_diagnosis_stably() {
    let baseline = mixed_preflight_failure(false);
    let reordered = mixed_preflight_failure(true);

    for failure in [&baseline, &reordered] {
        assert_eq!(
            failure.diagnosis(),
            ProblemDiagnosis::UninformativeSharedLevelSet,
            "unresolved semantics outrank a Direct Input Conflict"
        );
        assert!(failure.report().attempts().is_empty());
        assert_eq!(
            failure
                .report()
                .uninformative_shared_level_sets()
                .iter()
                .map(|evidence| (
                    evidence.group_id().as_str(),
                    evidence.member_source_id().as_str(),
                ))
                .collect::<Vec<_>>(),
            [("a-group", "z-member"), ("z-group", "a-member")]
        );
        let conflicts = failure.report().direct_input_conflicts();
        assert_eq!(conflicts.len(), 1);
        assert_eq!(conflicts[0].first_source().as_str(), "conflict-a");
        assert_eq!(conflicts[0].second_source().as_str(), "conflict-b");
        assert_eq!(conflicts[0].first_target(), 1.0);
        assert_eq!(conflicts[0].second_target(), 2.0);
    }

    assert_eq!(
        baseline.report().uninformative_shared_level_sets(),
        reordered.report().uninformative_shared_level_sets()
    );
    assert_eq!(
        baseline.report().direct_input_conflicts(),
        reordered.report().direct_input_conflicts()
    );
}

#[test]
fn cubic_pi1_rank_loss_is_tied_to_a_canonical_unidentified_field_mode() {
    let mut builder = problem_builder();
    builder
        .add(
            FieldValueObservation::try_new(SourceId::new("only-value"), point(1.0, 2.0, 3.0), 1.0)
                .unwrap(),
        )
        .unwrap();

    let failure = builder
        .build()
        .unwrap()
        .fit()
        .expect_err("one value does not identify the complete Cubic polynomial space");
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::UnidentifiedFieldMode);

    let interpreted = failure
        .report()
        .interpretable_rank_deficiency()
        .expect("the numerical rank loss is recovered to a canonical field concept");
    assert_eq!(
        interpreted.concept(),
        RankDeficiencyConcept::CubicPi1FieldMode
    );
    assert_eq!(
        interpreted.domain(),
        RankEvidenceDomain::CubicPolynomialPairing
    );
    assert_eq!(interpreted.source_ids(), &[SourceId::new("only-value")]);
    assert_eq!(
        interpreted
            .semantic_roles()
            .iter()
            .map(|role| role.as_str())
            .collect::<Vec<_>>(),
        ["field-value-observation/value"]
    );
    assert!(!interpreted.backend_invoked());
    assert!(!interpreted.hidden_regularization_applied());
    assert!(interpreted.canonical_mode_verified());
    assert!(interpreted.canonical_mode_residual() <= 1.0e-12);

    let rank = failure.report().rank_evidence().unwrap();
    assert_eq!(rank.domain(), interpreted.domain());
    assert_eq!(rank.decision(), RankDecision::RankDeficient);
}

#[test]
fn solve_attempt_termination_does_not_encode_canonical_acceptance() {
    let mut builder = problem_builder();
    for (source, location, value) in [
        ("origin", point(0.0, 0.0, 0.0), 1.0),
        ("east", point(1.0, 0.0, 0.0), 2.0),
        ("north", point(0.0, 1.0, 0.0), 3.0),
        ("up", point(0.0, 0.0, 1.0), 4.0),
    ] {
        builder
            .add(FieldValueObservation::try_new(SourceId::new(source), location, value).unwrap())
            .unwrap();
    }

    let success = builder.build().unwrap().fit().unwrap();
    let report = success.report();
    assert!(report.canonical_acceptance().unwrap().accepted());
    assert!(!report.attempts().is_empty());
    assert!(report.attempts().iter().all(|attempt| {
        attempt.termination() == SolveAttemptTermination::CandidateProduced
            && attempt.failure_reason().is_none()
    }));
}

#[test]
fn direct_conflict_outranks_gauge_evidence_without_discarding_it() {
    let mut builder = problem_builder();
    for (source, x_component) in [("gradient-a", 1.0), ("gradient-b", 2.0)] {
        builder
            .add(GradientObservation::new(
                SourceId::new(source),
                point(0.0, 0.0, 0.0),
                Vector3::try_new(x_component, 0.0, 0.0).unwrap(),
            ))
            .unwrap();
    }

    let failure = builder.build().unwrap().fit().unwrap_err();
    assert_eq!(failure.diagnosis(), ProblemDiagnosis::DirectInputConflict);
    assert!(failure.report().unidentified_additive_gauge().is_some());
    assert_eq!(failure.report().direct_input_conflicts().len(), 1);
    assert_eq!(
        failure.report().direct_input_conflicts()[0]
            .semantic_role()
            .as_str(),
        "gradient-observation/component/0"
    );
    assert!(failure.report().attempts().is_empty());
}

#[test]
fn unidentified_gauge_outranks_capacity_and_retains_the_checked_plan() {
    let mut builder = problem_builder();
    for index in 0..3_334 {
        builder
            .add(GradientObservation::new(
                SourceId::new(format!("gradient-{index:04}")),
                point(index as f64, 0.0, 0.0),
                Vector3::try_new(1.0, 0.0, 0.0).unwrap(),
            ))
            .unwrap();
    }

    let failure = builder.build().unwrap().fit().unwrap_err();
    assert_eq!(
        failure.diagnosis(),
        ProblemDiagnosis::UnidentifiedAdditiveGauge
    );
    let capacity = failure
        .report()
        .capacity()
        .expect("the lower-priority checked capacity evidence is retained");
    assert!(capacity.planned_peak_bytes().unwrap() > capacity.limit_bytes());
    assert!(!capacity.large_allocation_attempted());
    assert!(!capacity.backend_invocation_attempted());
    assert!(failure.report().attempts().is_empty());
}
