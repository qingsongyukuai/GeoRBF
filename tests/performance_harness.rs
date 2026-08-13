#[path = "../benches/support/mod.rs"]
mod support;

use georbf::{
    assemble_system, assemble_system_with_layout, constraint_layout, fit_single_surface_linear,
    Axis, DerivativePoint, FunctionalKernel, IsotropicKernel, ModelType, Parameters, Point,
    RbfKernel,
};
use support::{
    benchmark_case, evaluate_gradients_with_threads, evaluate_scalars_with_threads,
    result_checksum, BENCHMARK_CASE, FIXED_MULTI_THREADS,
};

#[test]
fn fixed_case_is_deterministic_and_release_sized() {
    let first = benchmark_case().expect("the fixed benchmark case is valid");
    let second = benchmark_case().expect("the fixed benchmark case is repeatable");

    assert_eq!(BENCHMARK_CASE, "single_surface_cubic_dense_v1");
    assert_eq!(FIXED_MULTI_THREADS, 2);
    assert_eq!(first.constraints.interfaces.len(), 96);
    assert_eq!(first.constraints.planars.len(), 16);
    assert_eq!(first.constraints.tangents.len(), 8);
    assert_eq!(first.queries.len(), 4_096);
    assert_eq!(first.checksum, second.checksum);
}

#[test]
fn fixed_thread_evaluation_preserves_order_and_bits() {
    let case = benchmark_case().expect("the fixed benchmark case is valid");
    let mut parameters = Parameters {
        model_type: ModelType::SingleSurface,
        ..Parameters::default()
    };
    parameters.polynomial_order = 1;
    let model = fit_single_surface_linear(&case.constraints, &parameters)
        .expect("the release benchmark case must fit");

    let scalar_single = evaluate_scalars_with_threads(&model, &case.queries, 1).unwrap();
    let scalar_multi =
        evaluate_scalars_with_threads(&model, &case.queries, FIXED_MULTI_THREADS).unwrap();
    let gradient_single = evaluate_gradients_with_threads(&model, &case.queries, 1).unwrap();
    let gradient_multi =
        evaluate_gradients_with_threads(&model, &case.queries, FIXED_MULTI_THREADS).unwrap();

    assert_eq!(scalar_single, scalar_multi);
    assert_eq!(gradient_single, gradient_multi);
    assert_eq!(
        result_checksum(&scalar_single, &gradient_single),
        result_checksum(&scalar_multi, &gradient_multi)
    );
}

#[test]
fn precomputed_layout_assembly_is_bit_identical() {
    let case = benchmark_case().expect("the fixed benchmark case is valid");
    let mut constraints = case.constraints.clone();
    constraints.remove_collocated();
    let kernel = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let layout = constraint_layout(ModelType::SingleSurface, &constraints, &case.parameters)
        .expect("the fixed layout is valid");
    let ordinary = assemble_system(
        &constraints,
        &case.parameters,
        FunctionalKernel::from(&kernel),
    )
    .unwrap();
    let precomputed = assemble_system_with_layout(
        layout,
        &constraints,
        &case.parameters,
        FunctionalKernel::from(&kernel),
    )
    .unwrap();

    assert_eq!(ordinary, precomputed);
}

#[test]
fn precomputed_layout_rejects_mismatched_constraint_counts() {
    let case = benchmark_case().expect("the fixed benchmark case is valid");
    let kernel = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let layout = constraint_layout(
        ModelType::SingleSurface,
        &case.constraints,
        &case.parameters,
    )
    .expect("the fixed layout is valid");
    let mut mismatched = case.constraints.clone();
    mismatched.interfaces.pop();

    assert!(assemble_system_with_layout(
        layout,
        &mismatched,
        &case.parameters,
        FunctionalKernel::from(&kernel),
    )
    .is_err());
}

#[test]
fn cached_cubic_spatial_derivatives_are_bit_identical() {
    let first = Point::new(1.25, -2.5, 3.75).unwrap();
    let second = Point::new(-4.5, 5.25, -6.0).unwrap();
    let kernel = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let cached_first = kernel
        .first_derivative_vector(&first, &second, DerivativePoint::First)
        .unwrap();
    let cached_second = kernel
        .first_derivative_vector(&first, &second, DerivativePoint::Second)
        .unwrap();
    let cached_hessian = kernel.mixed_hessian(&first, &second).unwrap();
    let axes = [Axis::X, Axis::Y, Axis::Z];

    for (index, axis) in axes.into_iter().enumerate() {
        assert_eq!(
            cached_first[index].to_bits(),
            kernel
                .first_derivative(&first, &second, DerivativePoint::First, axis)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            cached_second[index].to_bits(),
            kernel
                .first_derivative(&first, &second, DerivativePoint::Second, axis)
                .unwrap()
                .to_bits()
        );
        for (second_index, second_axis) in axes.into_iter().enumerate() {
            assert_eq!(
                cached_hessian[index][second_index].to_bits(),
                kernel
                    .mixed_second_derivative(&first, &second, axis, second_axis)
                    .unwrap()
                    .to_bits()
            );
        }
    }
}
