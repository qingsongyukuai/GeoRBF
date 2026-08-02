pub(crate) fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

pub(crate) fn dot3(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}
