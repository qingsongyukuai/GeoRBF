//! Surfe-compatible point ordering, collocation, and indexed value sorting.
//!
//! Sources:
//! - `surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`Point::operator<`, `collocated`)
//! - `surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`GRBF_Modelling_Methods::remove_collocated_constraints`)
//! - `math_lib/math_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//!   (`Math_methods::sort_vector_w_index`)

use std::cmp::Ordering;

use crate::{Constraints, Point, POSITION_EPSILON};

/// Compare only `(x, y, z)`, with frozen Surfe's lexicographic semantics.
///
/// GeoRBF constructors reject NaN and infinity, so this is a total ordering of
/// accepted positions. Signed zero compares equal, exactly as it does in the
/// C++ relational operators. A stable Rust sort supplies deterministic order
/// for points whose three coordinates compare equal.
pub fn compare_points(left: &Point, right: &Point) -> Ordering {
    for (left, right) in [
        (left.x(), right.x()),
        (left.y(), right.y()),
        (left.z(), right.z()),
    ] {
        if left < right {
            return Ordering::Less;
        }
        if left > right {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

/// Frozen Surfe's strict, axis-wise `1e-3` same-position predicate.
///
/// The fourth coordinate and all constraint payload fields are deliberately
/// ignored. Exactly `1e-3` on any axis is not collocated.
pub fn collocated(left: &Point, right: &Point) -> bool {
    (left.x() - right.x()).abs() < POSITION_EPSILON
        && (left.y() - right.y()).abs() < POSITION_EPSILON
        && (left.z() - right.z()).abs() < POSITION_EPSILON
}

/// Number removed from each independently cleaned constraint category.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CollocationRemoval {
    pub inequalities: usize,
    pub interfaces: usize,
    pub planars: usize,
    pub tangents: usize,
}

impl CollocationRemoval {
    pub const fn total(self) -> usize {
        self.inequalities + self.interfaces + self.planars + self.tangents
    }
}

impl Constraints {
    /// Sort each category by position and remove adjacent collocated values.
    ///
    /// Categories never remove values from one another. As in `std::unique`,
    /// the first value in each sorted collocation run supplies the retained
    /// level, normal, tangent, or other category-specific payload.
    pub fn remove_collocated(&mut self) -> CollocationRemoval {
        CollocationRemoval {
            inequalities: sort_and_remove(&mut self.inequalities, |value| value.point()),
            interfaces: sort_and_remove(&mut self.interfaces, |value| value.point()),
            planars: sort_and_remove(&mut self.planars, |value| value.point()),
            tangents: sort_and_remove(&mut self.tangents, |value| value.point()),
        }
    }
}

fn sort_and_remove<T>(values: &mut Vec<T>, point: impl Fn(&T) -> &Point + Copy) -> usize {
    values.sort_by(|left, right| compare_points(point(left), point(right)));
    let original_len = values.len();
    values.dedup_by(|left, right| collocated(point(left), point(right)));
    original_len - values.len()
}

/// Sort values ascending while applying the exact same permutation to indices.
///
/// This is a safe transliteration of frozen Surfe's Numerical Recipes
/// `sort_vector_w_index`, including its insertion-sort threshold, partition
/// scheduling, duplicate ordering, signed-zero ordering, and 50-slot stack.
/// It returns `false` without mutation for a length mismatch. Accepted GeoRBF
/// paths supply finite values, as required by the source algorithm.
pub fn sort_values_with_indices(values: &mut [f64], indices: &mut [usize]) -> bool {
    if values.len() != indices.len() {
        return false;
    }
    if values.is_empty() {
        return true;
    }

    const INSERTION_THRESHOLD: usize = 7;
    const STACK_SIZE: usize = 50;

    let mut stack = [0_usize; STACK_SIZE];
    let mut stack_top = -1_isize;
    let mut left = 0_usize;
    let mut right = values.len() - 1;

    loop {
        if right - left < INSERTION_THRESHOLD {
            for current in (left + 1)..=right {
                let value = values[current];
                let index = indices[current];
                let mut previous = current as isize - 1;
                while previous >= left as isize {
                    let previous_index = previous as usize;
                    if values[previous_index] <= value {
                        break;
                    }
                    values[previous_index + 1] = values[previous_index];
                    indices[previous_index + 1] = indices[previous_index];
                    previous -= 1;
                }
                let insertion = (previous + 1) as usize;
                values[insertion] = value;
                indices[insertion] = index;
            }

            if stack_top < 0 {
                break;
            }
            right = stack[stack_top as usize];
            stack_top -= 1;
            left = stack[stack_top as usize];
            stack_top -= 1;
        } else {
            let middle = (left + right) >> 1;
            values.swap(middle, left + 1);
            indices.swap(middle, left + 1);
            swap_if_greater(values, indices, left, right);
            swap_if_greater(values, indices, left + 1, right);
            swap_if_greater(values, indices, left, left + 1);

            let mut lower = left + 1;
            let mut upper = right;
            let pivot = values[left + 1];
            let pivot_index = indices[left + 1];
            loop {
                lower += 1;
                while values[lower] < pivot {
                    lower += 1;
                }
                upper -= 1;
                while values[upper] > pivot {
                    upper -= 1;
                }
                if upper < lower {
                    break;
                }
                values.swap(lower, upper);
                indices.swap(lower, upper);
            }

            values[left + 1] = values[upper];
            values[upper] = pivot;
            indices[left + 1] = indices[upper];
            indices[upper] = pivot_index;

            stack_top += 2;
            if stack_top as usize >= STACK_SIZE {
                return false;
            }
            if right - lower + 1 >= upper - left {
                stack[stack_top as usize] = right;
                stack[stack_top as usize - 1] = lower;
                right = upper - 1;
            } else {
                stack[stack_top as usize] = upper - 1;
                stack[stack_top as usize - 1] = left;
                left = lower;
            }
        }
    }

    true
}

fn swap_if_greater(values: &mut [f64], indices: &mut [usize], left: usize, right: usize) {
    if values[left] > values[right] {
        values.swap(left, right);
        indices.swap(left, right);
    }
}
