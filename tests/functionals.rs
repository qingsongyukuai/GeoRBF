use georbf::{
    AnisotropicKernel, Axis, DofLabel, FunctionalKernel, FunctionalPrimitive, Interface,
    IsotropicKernel, KernelError, LinearFunctional, ModifiedKernel, Planar, Point, RbfKernel,
    Tangent,
};

const AXES: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

fn point(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z).unwrap()
}

fn tangent(x: f64, y: f64, z: f64, tx: f64, ty: f64, tz: f64) -> Tangent {
    Tangent::new(x, y, z, tx, ty, tz).unwrap()
}

fn interface(x: f64, y: f64, z: f64) -> Interface {
    Interface::new(x, y, z, 20.0).unwrap()
}

fn groups() -> Vec<Vec<Interface>> {
    vec![vec![
        interface(0.0, 0.0, 0.0),
        interface(10.0, 1.0, 1.0),
        interface(2.0, -3.0, 0.0),
        interface(3.0, 4.0, 2.0),
        interface(4.0, 2.0, 0.5),
        interface(5.0, -1.0, 1.5),
    ]]
}

fn oblique_planars() -> Vec<Planar> {
    vec![
        Planar::from_normal(0.0, 0.0, 0.0, 1.0, 0.0, 0.0).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, 0.0, 1.0, 0.0).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, 0.36, -0.48, 0.8).unwrap(),
        Planar::from_normal(0.0, 0.0, 0.0, -0.48, 0.64, 0.6).unwrap(),
    ]
}

fn primitive_values(kernel: FunctionalKernel<'_>) -> Vec<f64> {
    let first = point(-1.25, 0.5, 2.0);
    let second = point(0.75, -0.25, 1.0);
    let first_tangent = tangent(-1.25, 0.5, 2.0, 0.3, -0.4, 0.5);
    let second_tangent = tangent(0.75, -0.25, 1.0, -0.2, 0.6, 0.7);
    let first_value = LinearFunctional::value(first.clone());
    let second_value = LinearFunctional::value(second.clone());
    let first_derivatives = AXES.map(|axis| LinearFunctional::derivative(first.clone(), axis));
    let second_derivatives = AXES.map(|axis| LinearFunctional::derivative(second.clone(), axis));
    let first_tangent = LinearFunctional::tangent(first_tangent);
    let second_tangent = LinearFunctional::tangent(second_tangent);
    let mut values = Vec::new();

    values.push(kernel.apply(&first_value, &second_value).unwrap());
    for derivative in &second_derivatives {
        values.push(kernel.apply(&first_value, derivative).unwrap());
    }
    for derivative in &first_derivatives {
        values.push(kernel.apply(derivative, &second_value).unwrap());
    }
    for first_derivative in &first_derivatives {
        for second_derivative in &second_derivatives {
            values.push(kernel.apply(first_derivative, second_derivative).unwrap());
        }
    }
    values.push(kernel.apply(&first_value, &second_tangent).unwrap());
    values.push(kernel.apply(&first_tangent, &second_value).unwrap());
    for derivative in &first_derivatives {
        values.push(kernel.apply(derivative, &second_tangent).unwrap());
    }
    for derivative in &second_derivatives {
        values.push(kernel.apply(&first_tangent, derivative).unwrap());
    }
    values.push(kernel.apply(&first_tangent, &second_tangent).unwrap());
    values
}

fn difference_values(kernel: FunctionalKernel<'_>) -> Vec<f64> {
    let row_positive = point(-0.5, 1.25, 0.75);
    let row_negative = point(1.5, -0.75, 2.25);
    let column_positive = point(0.25, -1.5, 1.75);
    let column_negative = point(2.0, 0.5, -0.25);
    let value_point = point(-1.0, 0.25, 1.5);
    let tangent = tangent(0.5, -0.25, 2.0, -0.35, 0.45, 0.8);
    let row_difference = LinearFunctional::difference(row_positive, row_negative);
    let column_difference = LinearFunctional::difference(column_positive, column_negative);
    let value = LinearFunctional::value(value_point.clone());
    let derivatives = AXES.map(|axis| LinearFunctional::derivative(value_point.clone(), axis));
    let tangent = LinearFunctional::tangent(tangent);
    let mut values = Vec::new();

    values.push(kernel.apply(&row_difference, &value).unwrap());
    values.push(kernel.apply(&value, &column_difference).unwrap());
    for derivative in &derivatives {
        values.push(kernel.apply(&row_difference, derivative).unwrap());
    }
    for derivative in &derivatives {
        values.push(kernel.apply(derivative, &column_difference).unwrap());
    }
    values.push(kernel.apply(&row_difference, &tangent).unwrap());
    values.push(kernel.apply(&tangent, &column_difference).unwrap());
    values.push(kernel.apply(&row_difference, &column_difference).unwrap());
    values
}

fn assert_bits(actual: &[f64], expected: &[u64]) {
    assert_eq!(actual.len(), expected.len());
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            actual.to_bits(),
            *expected,
            "value {index}: actual={actual:.17e}, expected bits={expected:016x}"
        );
    }
}

#[test]
fn every_primitive_combination_matches_frozen_cpp_golden() {
    const ISOTROPIC: [u64; 25] = [
        0x402a3cfe9be10415,
        0x402c4d4c2993cae2,
        0xc01539f91f2ed82a,
        0xc01c4d4c2993cae2,
        0xc02c4d4c2993cae2,
        0x401539f91f2ed82a,
        0x401c4d4c2993cae2,
        0xc02853b3155b1480,
        0x3ffe872701b38d2e,
        0x40045a1a01225e1e,
        0x3ffe872701b38d2e,
        0xc01f29f7d1bca01e,
        0xbfee872701b38d2e,
        0x40045a1a01225e1e,
        0xbfee872701b38d2e,
        0xc020b1e954ee3135,
        0xc025ef1b069f56d6,
        0xc006a43cee0fd582,
        0x40156ee362cb64b7,
        0xc016e55d4146a9e1,
        0xc01bb29f63f1927e,
        0xc0091f3819666c2d,
        0x4009b1c0d4a1ca06,
        0xc0083b46f62684dc,
        0x3fdbd4d05c6e5b98,
    ];
    const ANISOTROPIC: [u64; 25] = [
        0x403eb0cb0e38df50,
        0x403cedf1e0751366,
        0xc03491ea1f3fab2c,
        0xc032c90dd250b6c1,
        0xc03cedf1e0751366,
        0x403491ea1f3fab2c,
        0x4032c90dd250b6c1,
        0xc035dd258a3a2eaf,
        0x4025f181800960c2,
        0x40179c205e10b7d8,
        0x4025f181800960c2,
        0xc033ffe51594225e,
        0xc010c99bd6ba3394,
        0x40179c205e10b7d8,
        0xc010c99bd6ba3394,
        0xc0369e3e3556bf5b,
        0xc03f471352dc841c,
        0xc01e0e7d3450eeda,
        0x402e2cce71d5f379,
        0xc03121c4df5400ec,
        0xc03388048a3fe3ad,
        0xc01ffc8476de8fe8,
        0x402262c3f5312ba0,
        0xc01f70345ef79954,
        0x3ff9d050bae8ca70,
    ];
    const MODIFIED_ISOTROPIC: [u64; 25] = [
        0x4081366098816242,
        0xc069ff9a8c7a23dc,
        0xc062edf6108ea854,
        0x40851250e9a5ce9c,
        0xc04a4f64fe7ffdb2,
        0xc051fcbacd927580,
        0x406f71786eaa8fbb,
        0x404501b96ebe829e,
        0x40370d80bc4b03c7,
        0xc054bdd56fe55c80,
        0x4037a463df07288a,
        0x40390230709c7320,
        0xc054d1fbd4e42390,
        0xc054d25cc6f4b28a,
        0xc050d4b43b8f61f8,
        0x4072b3aead1aad7b,
        0x407a6bd08363258a,
        0x40615886cd5c2c57,
        0xc04a52dc7a50b416,
        0xc048028f382ac680,
        0x4067372e2792cdbc,
        0xc0433f85ebf02d10,
        0xc0425fea9c3d053a,
        0x4063c1343a1919eb,
        0x405811a9ba0f40d2,
    ];
    const MODIFIED_ANISOTROPIC: [u64; 25] = [
        0x408c4479384b9013,
        0xc0730fa15cb5b8bd,
        0xc06ca25a0e304e7b,
        0x4090f38406a20006,
        0xc050aec6a6372d8d,
        0xc06273be4d74f594,
        0x407a673b295e0326,
        0x404d7b67e0e199c4,
        0x4038bf84e5839e56,
        0xc05c73314a1b8690,
        0x403b88ae98335832,
        0x4045c15008b9881a,
        0xc0634ac5197e151a,
        0xc05e4c8e79ef04d1,
        0xc05b4a33ec003974,
        0x407e9c0fe3d20677,
        0x408557bb76d42056,
        0x406f481caf6aeb59,
        0xc05326a5c1894905,
        0xc055dbeb0f1505ec,
        0x4072d9240ef888f3,
        0xc04af6123b4f3520,
        0xc05023ccccee8c12,
        0x4071078bab0cf6f6,
        0x406458ed68282594,
    ];

    let isotropic = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let anisotropic = AnisotropicKernel::new(RbfKernel::Cubic, 1.0, &oblique_planars()).unwrap();
    let modified_isotropic = ModifiedKernel::from_isotropic(isotropic, &groups()).unwrap();
    let modified_anisotropic = ModifiedKernel::from_anisotropic(anisotropic, &groups()).unwrap();

    assert_bits(
        &primitive_values(FunctionalKernel::from(&isotropic)),
        &ISOTROPIC,
    );
    assert_bits(
        &primitive_values(FunctionalKernel::from(&anisotropic)),
        &ANISOTROPIC,
    );
    assert_bits(
        &primitive_values(FunctionalKernel::from(&modified_isotropic)),
        &MODIFIED_ISOTROPIC,
    );
    assert_bits(
        &primitive_values(FunctionalKernel::from(&modified_anisotropic)),
        &MODIFIED_ANISOTROPIC,
    );
}

#[test]
fn every_difference_combination_matches_increment_pair_golden() {
    const ISOTROPIC: [u64; 11] = [
        0xc033657da6c7dd78,
        0xc040092d321c353e,
        0x4032f1965ae63e47,
        0xc028d9275f9dc719,
        0x4022a2dd87b65553,
        0x4037384a69e93efe,
        0x402bf4d0e4979e77,
        0xc033e7a1729a5301,
        0xbff7410b6217fbba,
        0xc02c340e58f394fe,
        0x4038f5f8acf72c2f,
    ];
    const ANISOTROPIC: [u64; 11] = [
        0xc047217388b470bf,
        0xc04709520028b480,
        0x4045073d4804dcff,
        0xc0432665df88cf3d,
        0x4037ecf5b00892ce,
        0x403ab833664b8d99,
        0x4032da2cef5b25c2,
        0xc04869c5ee83f1da,
        0xc00f59afd9b1d128,
        0xc0423995aea73e24,
        0x4052e4764f2c8adc,
    ];
    const MODIFIED_ISOTROPIC: [u64; 11] = [
        0xc0867e14cd9fdca4,
        0x40980a9a9bbf80fe,
        0x405453b4b449ba1c,
        0x406066d5677dc33a,
        0xc07b74a6715ccdf0,
        0xc07167be5d9cb08a,
        0xc06e316f1de65881,
        0x408b3575b3b804ea,
        0xc074d983f54735c0,
        0x4085cb77814591be,
        0xc093c4b762e4886c,
    ];
    const MODIFIED_ANISOTROPIC: [u64; 11] = [
        0xc092eac3a5a1868c,
        0x40a2fa33e4794933,
        0x40566649b27f715c,
        0x40717a5497360ef9,
        0xc087fdfb851e1710,
        0xc07778896453a3bd,
        0xc079c6bcb56ea113,
        0x4095eb86e1af10d9,
        0xc0818bf5511911ba,
        0x40914dd367ed7eb9,
        0xc0a15251c3ef39e2,
    ];

    let isotropic = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let anisotropic = AnisotropicKernel::new(RbfKernel::Cubic, 1.0, &oblique_planars()).unwrap();
    let modified_isotropic = ModifiedKernel::from_isotropic(isotropic, &groups()).unwrap();
    let modified_anisotropic = ModifiedKernel::from_anisotropic(anisotropic, &groups()).unwrap();

    assert_bits(
        &difference_values(FunctionalKernel::from(&isotropic)),
        &ISOTROPIC,
    );
    assert_bits(
        &difference_values(FunctionalKernel::from(&anisotropic)),
        &ANISOTROPIC,
    );
    assert_bits(
        &difference_values(FunctionalKernel::from(&modified_isotropic)),
        &MODIFIED_ISOTROPIC,
    );
    assert_bits(
        &difference_values(FunctionalKernel::from(&modified_anisotropic)),
        &MODIFIED_ANISOTROPIC,
    );
}

#[test]
fn representations_have_unambiguous_labels_and_difference_expansion() {
    let location = point(1.0, 2.0, 3.0);
    let direction = tangent(4.0, 5.0, 6.0, -0.25, 0.5, 0.75);
    let value = LinearFunctional::value(location.clone());
    let derivative = LinearFunctional::derivative(location.clone(), Axis::Y);
    let tangent = LinearFunctional::tangent(direction);
    let difference = LinearFunctional::difference(location, point(-1.0, -2.0, -3.0));

    assert_eq!(value.label(), DofLabel::Value);
    assert_eq!(derivative.label(), DofLabel::Derivative(Axis::Y));
    assert_eq!(tangent.label(), DofLabel::Tangent);
    assert_eq!(difference.label(), DofLabel::Difference);
    assert_eq!(value.expansion().len(), 1);
    assert_eq!(derivative.expansion().len(), 1);
    assert_eq!(tangent.expansion().len(), 1);

    let terms = difference.expansion();
    assert_eq!(terms.len(), 2);
    assert_eq!(terms[0].coefficient(), 1.0);
    assert_eq!(terms[1].coefficient(), -1.0);
    match terms[0].primitive() {
        FunctionalPrimitive::Value(point) => assert_eq!(point.position(), [1.0, 2.0, 3.0]),
        _ => panic!("difference minuend must expand to Value"),
    }
    match terms[1].primitive() {
        FunctionalPrimitive::Value(point) => assert_eq!(point.position(), [-1.0, -2.0, -3.0]),
        _ => panic!("difference subtrahend must expand to Value"),
    }
}

#[test]
fn difference_is_linear_and_preserves_frozen_parenthesization() {
    let kernel = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let kernel = FunctionalKernel::from(&kernel);
    let positive = point(-0.5, 1.25, 0.75);
    let negative = point(1.5, -0.75, 2.25);
    let column_positive = point(0.25, -1.5, 1.75);
    let column_negative = point(2.0, 0.5, -0.25);
    let row_difference = LinearFunctional::difference(positive.clone(), negative.clone());
    let column_difference =
        LinearFunctional::difference(column_positive.clone(), column_negative.clone());
    let row_positive = LinearFunctional::value(positive);
    let row_negative = LinearFunctional::value(negative);
    let column_positive = LinearFunctional::value(column_positive);
    let column_negative = LinearFunctional::value(column_negative);

    let actual = kernel.apply(&row_difference, &column_difference).unwrap();
    let v1 = kernel.apply(&row_positive, &column_positive).unwrap();
    let v2 = kernel.apply(&row_positive, &column_negative).unwrap();
    let v3 = kernel.apply(&row_negative, &column_positive).unwrap();
    let v4 = kernel.apply(&row_negative, &column_negative).unwrap();
    assert_eq!(actual.to_bits(), ((v1 - v2) - (v3 - v4)).to_bits());
}

#[test]
fn parameter_swap_and_radial_sign_identities_hold() {
    let kernel = IsotropicKernel::new(RbfKernel::Gaussian, 0.75);
    let kernel = FunctionalKernel::from(&kernel);
    let first = point(-0.75, 1.25, 0.5);
    let second = point(1.5, -0.5, 2.0);
    let value_first = LinearFunctional::value(first.clone());
    let value_second = LinearFunctional::value(second.clone());

    for axis in AXES {
        let derivative_first = LinearFunctional::derivative(first.clone(), axis);
        let derivative_second = LinearFunctional::derivative(second.clone(), axis);
        let left = kernel.apply(&value_first, &derivative_second).unwrap();
        let swapped = kernel.apply(&derivative_second, &value_first).unwrap();
        let opposite = kernel.apply(&derivative_first, &value_second).unwrap();
        assert_eq!(left.to_bits(), swapped.to_bits());
        assert_eq!(left.to_bits(), (-opposite).to_bits());
    }

    let first_tangent = LinearFunctional::tangent(tangent(-0.75, 1.25, 0.5, 0.2, -0.4, 0.7));
    let second_tangent = LinearFunctional::tangent(tangent(1.5, -0.5, 2.0, -0.3, 0.6, 0.1));
    let forward = kernel.apply(&first_tangent, &second_tangent).unwrap();
    let swapped = kernel.apply(&second_tangent, &first_tangent).unwrap();
    assert_eq!(forward.to_bits(), swapped.to_bits());
}

fn displaced(point: &Point, direction: [f64; 3], amount: f64) -> Point {
    Point::with_c(
        point.x() + direction[0] * amount,
        point.y() + direction[1] * amount,
        point.z() + direction[2] * amount,
        point.c(),
    )
    .unwrap()
}

fn assert_close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    let delta = (actual - expected).abs();
    let limit = absolute + relative * actual.abs().max(expected.abs());
    assert!(
        delta <= limit,
        "actual={actual:.17e}, expected={expected:.17e}, delta={delta:.3e}, limit={limit:.3e}"
    );
}

#[test]
fn derivative_and_tangent_functionals_match_finite_differences() {
    let radial = IsotropicKernel::new(RbfKernel::Gaussian, 0.6);
    let kernel = FunctionalKernel::from(&radial);
    let first = point(-0.75, 0.25, 1.5);
    let second = point(1.25, -0.5, 0.75);
    let first_value = LinearFunctional::value(first.clone());
    let second_value = LinearFunctional::value(second.clone());
    let step = 1.0e-5;

    for axis in AXES {
        let direction = match axis {
            Axis::X => [1.0, 0.0, 0.0],
            Axis::Y => [0.0, 1.0, 0.0],
            Axis::Z => [0.0, 0.0, 1.0],
        };
        let plus = LinearFunctional::value(displaced(&first, direction, step));
        let minus = LinearFunctional::value(displaced(&first, direction, -step));
        let finite = (kernel.apply(&plus, &second_value).unwrap()
            - kernel.apply(&minus, &second_value).unwrap())
            / (2.0 * step);
        let analytic = kernel
            .apply(
                &LinearFunctional::derivative(first.clone(), axis),
                &second_value,
            )
            .unwrap();
        assert_close(analytic, finite, 2.0e-10, 2.0e-9);

        let plus = LinearFunctional::value(displaced(&second, direction, step));
        let minus = LinearFunctional::value(displaced(&second, direction, -step));
        let finite = (kernel.apply(&first_value, &plus).unwrap()
            - kernel.apply(&first_value, &minus).unwrap())
            / (2.0 * step);
        let analytic = kernel
            .apply(
                &first_value,
                &LinearFunctional::derivative(second.clone(), axis),
            )
            .unwrap();
        assert_close(analytic, finite, 2.0e-10, 2.0e-9);
    }

    let direction = [0.3, -0.4, 0.5];
    let tangent = LinearFunctional::tangent(
        Tangent::new(
            first.x(),
            first.y(),
            first.z(),
            direction[0],
            direction[1],
            direction[2],
        )
        .unwrap(),
    );
    let plus = LinearFunctional::value(displaced(&first, direction, step));
    let minus = LinearFunctional::value(displaced(&first, direction, -step));
    let finite = (kernel.apply(&plus, &second_value).unwrap()
        - kernel.apply(&minus, &second_value).unwrap())
        / (2.0 * step);
    let analytic = kernel.apply(&tangent, &second_value).unwrap();
    assert_close(analytic, finite, 2.0e-10, 2.0e-9);

    let component_sum = direction[0]
        * kernel
            .apply(
                &LinearFunctional::derivative(first.clone(), Axis::X),
                &second_value,
            )
            .unwrap()
        + direction[1]
            * kernel
                .apply(
                    &LinearFunctional::derivative(first.clone(), Axis::Y),
                    &second_value,
                )
                .unwrap()
        + direction[2]
            * kernel
                .apply(&LinearFunctional::derivative(first, Axis::Z), &second_value)
                .unwrap();
    assert_eq!(analytic.to_bits(), component_sum.to_bits());

    let second_direction = [-0.2, 0.6, 0.7];
    let second_tangent = LinearFunctional::tangent(
        Tangent::new(
            second.x(),
            second.y(),
            second.z(),
            second_direction[0],
            second_direction[1],
            second_direction[2],
        )
        .unwrap(),
    );
    let plus = LinearFunctional::value(displaced(&second, second_direction, step));
    let minus = LinearFunctional::value(displaced(&second, second_direction, -step));
    let finite = (kernel.apply(&first_value, &plus).unwrap()
        - kernel.apply(&first_value, &minus).unwrap())
        / (2.0 * step);
    let analytic = kernel.apply(&first_value, &second_tangent).unwrap();
    assert_close(analytic, finite, 2.0e-10, 2.0e-9);
}

#[test]
fn linear_radius_errors_propagate_through_all_composites() {
    let radial = IsotropicKernel::new(RbfKernel::Linear, 1.0);
    let kernel = FunctionalKernel::from(&radial);
    let first = point(0.0, 0.0, 0.0);
    let second = point(1.0, 2.0, 3.0);
    let value = LinearFunctional::value(first.clone());
    let derivative = LinearFunctional::derivative(second.clone(), Axis::X);
    let difference = LinearFunctional::difference(first, second);

    assert_eq!(
        kernel.apply(&value, &derivative),
        Err(KernelError::LinearDerivativeUnavailable)
    );
    assert_eq!(
        kernel.apply(&difference, &derivative),
        Err(KernelError::LinearDerivativeUnavailable)
    );
}
