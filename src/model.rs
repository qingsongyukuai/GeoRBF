//! Immutable solved models and single-point field sampling.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;
use std::sync::Arc;

use crate::cubic_equality::RecoveredCubicField;
use crate::functional::GroupId;
use crate::geometry::{Point3, Vector3};
use crate::problem::ProblemSnapshot;

/// An owning, immutable, cheaply cloned solved field model.
#[derive(Clone)]
pub struct SolvedModel {
    inner: Arc<SolvedModelData>,
}

impl fmt::Debug for SolvedModel {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SolvedModel { .. }")
    }
}

impl SolvedModel {
    pub(crate) fn new(
        snapshot: ProblemSnapshot,
        field: RecoveredCubicField,
        shared_level_values: BTreeMap<GroupId, f64>,
    ) -> Self {
        Self {
            inner: Arc::new(SolvedModelData {
                snapshot,
                field,
                shared_level_values,
            }),
        }
    }

    /// Evaluates field value and complete gradient together at one point.
    pub fn evaluate(&self, point: Point3) -> Result<FieldSample, QueryError> {
        let sample = self.inner.field.sample(point.components());
        if !sample.value.is_finite() || sample.gradient.iter().any(|value| !value.is_finite()) {
            return Err(QueryError {
                reason: QueryErrorReason::NonFiniteResult,
            });
        }
        let gradient = Vector3::try_new(sample.gradient[0], sample.gradient[1], sample.gradient[2])
            .expect("the query finiteness check covers every gradient component");
        Ok(FieldSample {
            value: sample.value,
            gradient,
        })
    }

    /// Returns the immutable problem snapshot this model permanently solves.
    pub fn problem_snapshot(&self) -> &ProblemSnapshot {
        &self.inner.snapshot
    }

    /// Returns a recovered shared-level or horizon value by stable GroupId.
    pub fn shared_level_value(&self, group_id: &GroupId) -> Option<f64> {
        self.inner.shared_level_values.get(group_id).copied()
    }
}

#[derive(Debug)]
struct SolvedModelData {
    snapshot: ProblemSnapshot,
    field: RecoveredCubicField,
    shared_level_values: BTreeMap<GroupId, f64>,
}

/// One coherent field value and complete input-frame gradient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FieldSample {
    value: f64,
    gradient: Vector3,
}

impl FieldSample {
    /// Returns the scalar field value.
    pub fn value(self) -> f64 {
        self.value
    }

    /// Returns the complete gradient in the declared input coordinate frame.
    pub fn gradient(self) -> Vector3 {
        self.gradient
    }
}

/// Why a model query was rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum QueryErrorReason {
    /// Evaluation overflowed or otherwise produced a non-finite observable.
    NonFiniteResult,
}

/// Atomic failure of a model query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QueryError {
    reason: QueryErrorReason,
}

impl QueryError {
    /// Returns the structured query-failure reason.
    pub fn reason(self) -> QueryErrorReason {
        self.reason
    }
}

impl fmt::Display for QueryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("model evaluation did not produce finite observables")
    }
}

impl Error for QueryError {}
