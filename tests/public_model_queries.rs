use std::num::NonZeroUsize;
use std::sync::Arc;
use std::thread;

use georbf::fit::FitSuccess;
use georbf::geometry::{
    FieldUnitLabel, Handedness, InputCoordinateFrame, LengthUnitLabel, Point3, Vector3,
};
use georbf::kernel::KernelKind;
use georbf::model::{FieldSample, QueryErrorReason};
use georbf::observation::{FieldValueObservation, GradientObservation};
use georbf::problem::{FitConfiguration, ThreadBudget};
use georbf::{ProblemBuilder, ProblemSnapshot, SolvedModel, SourceId};

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::try_new(x, y, z).expect("the manufactured query point is finite")
}

fn affine_value(location: Point3) -> f64 {
    let [x, y, z] = location.components();
    2.0 + 0.5 * x - 1.25 * y + 0.75 * z
}

fn build_affine_snapshot(thread_budget: ThreadBudget) -> ProblemSnapshot {
    let frame = InputCoordinateFrame::try_new(
        ["east", "north", "elevation"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )
    .unwrap();
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("field-unit"));
    builder.set_fit_configuration(FitConfiguration::default().with_thread_budget(thread_budget));
    for (index, location) in [
        point(-1.0, -1.0, -1.0),
        point(1.0, -1.0, -1.0),
        point(-1.0, 1.0, -1.0),
        point(-1.0, -1.0, 1.0),
        point(1.0, 1.0, 0.5),
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("value-{index}")),
                    location,
                    affine_value(location),
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder
        .add(GradientObservation::new(
            SourceId::new("gradient"),
            point(0.25, -0.5, 0.75),
            Vector3::try_new(0.5, -1.25, 0.75).unwrap(),
        ))
        .unwrap();
    builder.build().unwrap()
}

fn fit_affine_model() -> SolvedModel {
    build_affine_snapshot(ThreadBudget::Automatic)
        .fit()
        .unwrap()
        .into_parts()
        .0
}

fn fit_cancellation_model() -> FitSuccess {
    let frame = InputCoordinateFrame::try_new(
        ["x", "y", "z"],
        Handedness::Right,
        LengthUnitLabel::new("m"),
    )
    .unwrap();
    let mut builder = ProblemBuilder::new(frame, FieldUnitLabel::new("field-unit"));
    for (index, (location, value)) in [
        (point(0.0, 0.0, 0.0), 2.196_152_422_706_632),
        (point(1.0, 0.0, 0.0), -0.828_427_124_746_190_1),
        (point(0.0, 1.0, 0.0), -0.828_427_124_746_190_1),
        (point(0.0, 0.0, 1.0), -0.828_427_124_746_190_1),
        (point(1.0, 1.0, 1.0), 1.907_023_471_174_693_5),
    ]
    .into_iter()
    .enumerate()
    {
        builder
            .add(
                FieldValueObservation::try_new(
                    SourceId::new(format!("cancellation-{index}")),
                    location,
                    value,
                )
                .unwrap(),
            )
            .unwrap();
    }
    builder.build().unwrap().fit().unwrap()
}

fn assert_samples_equivalent(actual: FieldSample, expected: FieldSample) {
    let sample_reference_scale = actual.value().abs().max(expected.value().abs()).max(
        actual
            .gradient()
            .components()
            .into_iter()
            .chain(expected.gradient().components())
            .map(f64::abs)
            .fold(0.0, f64::max),
    );
    let field_scale = actual.value().abs().max(expected.value().abs()).max(1.0);
    let tolerance = 1.0e-12 * field_scale + 1.0e-11 * sample_reference_scale;
    assert!(
        (actual.value() - expected.value()).abs() <= tolerance,
        "value actual={:e}, expected={:e}, tolerance={tolerance:e}",
        actual.value(),
        expected.value()
    );
    for (actual, expected) in actual
        .gradient()
        .components()
        .into_iter()
        .zip(expected.gradient().components())
    {
        assert!(
            (actual - expected).abs() <= tolerance,
            "gradient actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
        );
    }
}

#[test]
fn logical_batch_matches_single_queries_in_input_order() {
    let model = fit_affine_model();
    let queries = [
        point(0.2, -0.3, 0.4),
        point(-3.0, 2.0, 0.5),
        point(1.5, -2.5, 4.0),
        point(0.0, 0.0, 0.0),
    ];
    let expected = queries
        .iter()
        .copied()
        .map(|query| model.evaluate(query).unwrap())
        .collect::<Vec<_>>();

    let actual = model
        .evaluate_batch(&queries)
        .expect("a finite logical query batch succeeds atomically");

    assert_eq!(actual, expected);
}

#[test]
fn logical_batch_failure_is_atomic_and_identifies_the_first_invalid_index() {
    let model = fit_affine_model();
    let queries = [
        point(0.2, -0.3, 0.4),
        point(f64::MAX, 0.0, 0.0),
        point(-f64::MAX, 0.0, 0.0),
    ];

    let failure = model
        .evaluate_batch(&queries)
        .expect_err("a non-finite field observable rejects the whole logical batch");

    assert_eq!(failure.reason(), QueryErrorReason::NonFiniteResult);
    assert_eq!(failure.point_index(), Some(1));
}

#[test]
fn empty_repeated_and_large_finite_queries_have_deterministic_batch_semantics() {
    let model = fit_affine_model();
    assert!(model.evaluate_batch(&[]).unwrap().is_empty());

    let repeated = point(0.25, -0.5, 0.75);
    let large = point(1.0e40, -1.0e40, 5.0e39);
    let queries = [repeated, repeated, large, repeated];
    let batch = model
        .evaluate_batch(&queries)
        .expect("large finite coordinates with finite observables remain valid queries");

    assert_eq!(batch[0], batch[1]);
    assert_eq!(batch[0], batch[3]);
    for (sample, query) in batch.into_iter().zip(queries) {
        assert_eq!(sample, model.evaluate(query).unwrap());
    }
}

#[test]
fn logical_batch_smoke_accepts_one_hundred_thousand_locations() {
    let model = fit_affine_model();
    let queries = (0..100_000)
        .map(|index| {
            let x = (index % 101) as f64 * 0.01 - 0.5;
            let y = (index % 37) as f64 * -0.02 + 0.25;
            let z = (index % 19) as f64 * 0.03 - 0.2;
            point(x, y, z)
        })
        .collect::<Vec<_>>();

    let samples = model
        .evaluate_batch(&queries)
        .expect("the guaranteed logical batch fits checked query scratch");

    assert_eq!(samples.len(), 100_000);
    for index in [0, 2_047, 2_048, 99_999] {
        assert_eq!(samples[index], model.evaluate(queries[index]).unwrap());
    }
}

#[test]
fn cloned_models_support_concurrent_single_and_varied_batch_queries() {
    fn assert_send_sync_clone<T: Send + Sync + Clone>() {}
    assert_send_sync_clone::<SolvedModel>();

    let model = fit_affine_model();
    let queries = Arc::new(
        (0..4_103)
            .map(|index| {
                let x = (index % 53) as f64 * 0.02 - 0.4;
                let y = (index % 29) as f64 * -0.03 + 0.3;
                let z = (index % 17) as f64 * 0.01 - 0.1;
                point(x, y, z)
            })
            .collect::<Vec<_>>(),
    );
    let expected = model.evaluate_batch(queries.as_slice()).unwrap();

    let handles = [1, 2_047, 2_049, 4_103].map(|length| {
        let model = model.clone();
        let queries = Arc::clone(&queries);
        thread::spawn(move || {
            let batch = model.evaluate_batch(&queries[..length]).unwrap();
            let single = model.evaluate(queries[length - 1]).unwrap();
            (batch, single)
        })
    });

    for (length, handle) in [1, 2_047, 2_049, 4_103].into_iter().zip(handles) {
        let (batch, single) = handle
            .join()
            .expect("read-only model query thread succeeds");
        assert_eq!(batch, expected[..length]);
        assert_eq!(single, expected[length - 1]);
    }
}

#[test]
fn model_retains_its_problem_contract_across_query_and_fit_resource_plans() {
    let automatic_snapshot = build_affine_snapshot(ThreadBudget::Automatic);
    let exact_snapshot = build_affine_snapshot(ThreadBudget::Exact(NonZeroUsize::new(1).unwrap()));
    let automatic = automatic_snapshot.fit().unwrap();
    let exact = exact_snapshot.fit().unwrap();
    let model = automatic.model();
    let owned_snapshot = model.problem_snapshot();

    assert_eq!(owned_snapshot.resolved_kernel().kind(), KernelKind::Cubic);
    assert_eq!(
        owned_snapshot
            .fit_configuration()
            .numerical_policy()
            .as_str(),
        "georbf-v2"
    );
    assert_eq!(
        owned_snapshot.input_coordinate_frame().axis_labels(),
        ["east", "north", "elevation"]
    );
    assert_eq!(
        owned_snapshot
            .input_coordinate_frame()
            .length_unit()
            .as_str(),
        "m"
    );
    assert_eq!(owned_snapshot.field_unit().as_str(), "field-unit");

    let queries = [point(0.2, -0.3, 0.4), point(-0.5, 0.25, 0.75)];
    let before = model.evaluate_batch(&queries).unwrap();
    let exact_samples = exact.model().evaluate_batch(&queries).unwrap();
    let after = model.evaluate_batch(&queries).unwrap();
    assert_eq!(before, after);
    for (actual, expected) in exact_samples.into_iter().zip(before) {
        assert_samples_equivalent(actual, expected);
    }
    assert_eq!(
        automatic.report().problem_size().center_coefficients(),
        exact.report().problem_size().center_coefficients()
    );
}

#[test]
fn recovered_field_publishes_verified_query_representation_evidence() {
    let fit = fit_cancellation_model();
    let representation = fit.report().representation_evidence();
    let query = representation
        .verified_query_representation()
        .expect("a successful fit publishes its verified query representation");

    assert!(query.verified());
    assert!(query.all_source_response_verified());
    assert!(query.pi1_side_condition_verified());
    assert!(query.field_energy_round_trip_verified());
    assert!(query.basis_round_trip_error() <= 1.0e-11);
}

#[test]
fn cancellation_queries_escalate_the_same_field_and_remain_batch_equivalent() {
    let fit = fit_cancellation_model();
    let model = fit.model();
    let query = point(100.0, 100.0, 100.0);
    let expected_value: f64 = 516.725_594_529_481_5;
    let expected_gradient: [f64; 3] = [1.732_060_494_437_222_7; 3];

    let single = model
        .evaluate(query)
        .expect("precision escalation reliably sums the verified field");
    let batch = model.evaluate_batch(&[query, query]).unwrap();

    assert_eq!(batch, vec![single, single]);
    let value_tolerance = 1.0e-12 * expected_value.abs().max(1.0) + 1.0e-11 * expected_value.abs();
    assert!((single.value() - expected_value).abs() <= value_tolerance);
    for (component, expected) in single
        .gradient()
        .components()
        .into_iter()
        .zip(expected_gradient)
    {
        let tolerance = 1.0e-12 * expected_value.abs().max(1.0) + 1.0e-11 * expected.abs();
        assert!((component - expected).abs() <= tolerance);
    }
}

#[test]
fn unresolved_query_cancellation_is_structured_and_batch_atomic() {
    let model = fit_cancellation_model().into_parts().0;
    let safe = point(0.25, 0.25, 0.25);
    let unresolved = point(1.0e110, 1.0e110, 1.0e110);

    let single = model
        .evaluate(unresolved)
        .expect_err("bounded precision cannot certify the cancelled observables");
    assert_eq!(single.reason(), QueryErrorReason::NumericalIndeterminate);
    assert_eq!(single.point_index(), None);

    let batch = model
        .evaluate_batch(&[safe, unresolved, safe])
        .expect_err("one indeterminate point rejects the logical batch atomically");
    assert_eq!(batch.reason(), QueryErrorReason::NumericalIndeterminate);
    assert_eq!(batch.point_index(), Some(1));
}
