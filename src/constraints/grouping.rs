//! Exact Surfe level and interface-reference grouping.
//!
//! Source: `surfe_lib/modeling_methods.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`_get_distinct_interface_iso_values`, `_get_interface_points`,
//! `_get_distinct_inequality_iso_values`, `get_interface_data`).

use std::cmp::Ordering;

use super::Constraints;

/// The three interface containers populated by frozen `get_interface_data`.
#[derive(Clone, Debug, PartialEq)]
pub struct InterfaceGrouping {
    levels_descending: Vec<f64>,
    reference_indices: Vec<usize>,
    multi_point_groups: Vec<Vec<usize>>,
}

impl InterfaceGrouping {
    /// Every exact interface level, largest to smallest.
    pub fn levels_descending(&self) -> &[f64] {
        &self.levels_descending
    }

    /// First point at each level, aligned with [`Self::levels_descending`].
    pub fn reference_indices(&self) -> &[usize] {
        &self.reference_indices
    }

    /// Per-level point indices for levels containing at least two points.
    ///
    /// Frozen Surfe drops singleton lists after recording their references.
    pub fn multi_point_groups(&self) -> &[Vec<usize>] {
        &self.multi_point_groups
    }

    /// Same-level difference degrees of freedom available to increment models.
    pub fn increment_pair_count(&self) -> usize {
        self.multi_point_groups
            .iter()
            .map(|group| group.len() - 1)
            .sum()
    }

    /// Adjacent reference pairs available to stratigraphic sequencing.
    pub fn sequenced_reference_pair_count(&self) -> usize {
        self.reference_indices.len().saturating_sub(1)
    }
}

impl Constraints {
    /// Group the current interface order by exact level.
    ///
    /// The public Surfe pipeline calls collocation removal first; callers that
    /// need that pipeline order should call [`Constraints::remove_collocated`]
    /// before this method. Empty interface input corresponds to the source
    /// `false` return and is represented by `None`.
    pub fn interface_grouping(&self) -> Option<InterfaceGrouping> {
        if self.interfaces.is_empty() {
            return None;
        }

        let levels_descending =
            distinct_levels_descending(self.interfaces.iter().map(|interface| interface.level()));
        let mut all_groups = Vec::with_capacity(levels_descending.len());
        for level in &levels_descending {
            all_groups.push(
                self.interfaces
                    .iter()
                    .enumerate()
                    .filter_map(|(index, interface)| (interface.level() == *level).then_some(index))
                    .collect::<Vec<_>>(),
            );
        }

        let reference_indices = all_groups.iter().map(|group| group[0]).collect();
        let multi_point_groups = all_groups
            .into_iter()
            .filter(|group| group.len() > 1)
            .collect();
        Some(InterfaceGrouping {
            levels_descending,
            reference_indices,
            multi_point_groups,
        })
    }

    /// Exact distinct inequality levels, largest to smallest.
    pub fn distinct_inequality_levels(&self) -> Vec<f64> {
        distinct_levels_descending(
            self.inequalities
                .iter()
                .map(|inequality| inequality.level()),
        )
    }
}

fn distinct_levels_descending(levels: impl IntoIterator<Item = f64>) -> Vec<f64> {
    let mut distinct = Vec::new();
    for level in levels {
        if !distinct.iter().any(|existing| *existing == level) {
            distinct.push(level);
        }
    }
    distinct.sort_by(|left, right| compare_finite(*right, *left));
    distinct
}

fn compare_finite(left: f64, right: f64) -> Ordering {
    if left < right {
        Ordering::Less
    } else if left > right {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}
