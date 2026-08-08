//! Immutable solved models and atomic single/batch field sampling.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;
use std::mem::size_of;
use std::sync::Arc;

use crate::cubic_equality::{
    FieldSample as RecoveredFieldSample, QuerySampleFailure, RecoveredCubicField,
};
use crate::functional::GroupId;
use crate::geometry::{Point3, Vector3};
use crate::problem::ProblemSnapshot;

const QUERY_SCRATCH_LIMIT_BYTES: u64 = 256 * 1024 * 1024;
const QUERY_CHUNK_TARGET_BYTES: u64 = 64 * 1024;
type RecoveredQueryResult = Result<RecoveredFieldSample, QuerySampleFailure>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct QueryPlan {
    output_capacity: usize,
    scratch_bytes: u64,
    chunk_len: usize,
    chunk_scratch_bytes: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QueryPlanningFailure {
    ArithmeticOverflow,
    ScratchLimitExceeded {
        planned_scratch_bytes: u64,
        limit_bytes: u64,
    },
}

fn plan_query_capacity(point_count: usize) -> Result<QueryPlan, QueryPlanningFailure> {
    let point_count =
        u64::try_from(point_count).map_err(|_| QueryPlanningFailure::ArithmeticOverflow)?;
    let output_bytes = point_count
        .checked_mul(size_of::<FieldSample>() as u64)
        .ok_or(QueryPlanningFailure::ArithmeticOverflow)?;
    if point_count == 0 {
        return Ok(QueryPlan {
            output_capacity: 0,
            scratch_bytes: 0,
            chunk_len: 0,
            chunk_scratch_bytes: 0,
        });
    }

    let bytes_per_recovered_sample = size_of::<RecoveredQueryResult>() as u64;
    let minimum_scratch_bytes = output_bytes
        .checked_add(bytes_per_recovered_sample)
        .ok_or(QueryPlanningFailure::ArithmeticOverflow)?;
    if minimum_scratch_bytes > QUERY_SCRATCH_LIMIT_BYTES {
        return Err(QueryPlanningFailure::ScratchLimitExceeded {
            planned_scratch_bytes: minimum_scratch_bytes,
            limit_bytes: QUERY_SCRATCH_LIMIT_BYTES,
        });
    }

    let available_chunk_bytes = QUERY_SCRATCH_LIMIT_BYTES - output_bytes;
    let maximum_chunk_len = available_chunk_bytes / bytes_per_recovered_sample;
    let target_chunk_len = (QUERY_CHUNK_TARGET_BYTES / bytes_per_recovered_sample)
        .max(1)
        .min(maximum_chunk_len);
    let target_chunk_len = usize::try_from(target_chunk_len)
        .map_err(|_| QueryPlanningFailure::ArithmeticOverflow)?
        .max(1);
    let point_count =
        usize::try_from(point_count).map_err(|_| QueryPlanningFailure::ArithmeticOverflow)?;
    let chunk_len = point_count.min(target_chunk_len);
    let chunk_scratch_bytes = u64::try_from(chunk_len)
        .map_err(|_| QueryPlanningFailure::ArithmeticOverflow)?
        .checked_mul(bytes_per_recovered_sample)
        .ok_or(QueryPlanningFailure::ArithmeticOverflow)?;
    let scratch_bytes = output_bytes
        .checked_add(chunk_scratch_bytes)
        .ok_or(QueryPlanningFailure::ArithmeticOverflow)?;
    Ok(QueryPlan {
        output_capacity: point_count,
        scratch_bytes,
        chunk_len,
        chunk_scratch_bytes,
    })
}

fn checked_field_sample(sample: RecoveredFieldSample) -> Result<FieldSample, QueryError> {
    if !sample.value.is_finite() || sample.gradient.iter().any(|value| !value.is_finite()) {
        return Err(QueryError {
            reason: QueryErrorReason::NonFiniteResult,
            point_index: None,
            planned_scratch_bytes: None,
            scratch_limit_bytes: None,
        });
    }
    let gradient = Vector3::try_new(sample.gradient[0], sample.gradient[1], sample.gradient[2])
        .expect("the query finiteness check covers every gradient component");
    Ok(FieldSample {
        value: sample.value,
        gradient,
    })
}

fn query_field_sample(
    field: &RecoveredCubicField,
    point: Point3,
) -> Result<FieldSample, QueryError> {
    let sample = field
        .reliable_sample(point.components())
        .map_err(QueryError::from_sample_failure)?;
    checked_field_sample(sample)
}

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
        query_field_sample(&self.inner.field, point)
    }

    /// Evaluates an ordered logical batch of points atomically.
    pub fn evaluate_batch(&self, points: &[Point3]) -> Result<Vec<FieldSample>, QueryError> {
        let plan = plan_query_capacity(points.len()).map_err(QueryError::from_planning_failure)?;
        debug_assert!(plan.scratch_bytes <= QUERY_SCRATCH_LIMIT_BYTES);
        let mut samples = Vec::with_capacity(plan.output_capacity);
        if plan.chunk_len == 0 {
            return Ok(samples);
        }
        for chunk in points.chunks(plan.chunk_len) {
            let recovered_samples = chunk
                .iter()
                .map(|point| self.inner.field.reliable_sample(point.components()))
                .collect::<Vec<_>>();
            debug_assert!(
                (recovered_samples.len() * size_of::<RecoveredQueryResult>()) as u64
                    <= plan.chunk_scratch_bytes
            );
            for sample in recovered_samples {
                let index = samples.len();
                let sample = sample
                    .map_err(QueryError::from_sample_failure)
                    .and_then(checked_field_sample)
                    .map_err(|error| error.with_point_index(index))?;
                samples.push(sample);
            }
        }
        Ok(samples)
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
    /// Bounded precision escalation could not certify a finite observable.
    NumericalIndeterminate,
    /// Checked atomic-batch scratch planning exceeded its resource envelope.
    CapacityExceeded,
}

/// Atomic failure of a model query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QueryError {
    reason: QueryErrorReason,
    point_index: Option<usize>,
    planned_scratch_bytes: Option<u64>,
    scratch_limit_bytes: Option<u64>,
}

impl QueryError {
    /// Returns the structured query-failure reason.
    pub fn reason(self) -> QueryErrorReason {
        self.reason
    }

    /// Returns the first failing logical-batch index, or `None` for a
    /// single-point query.
    pub fn point_index(self) -> Option<usize> {
        self.point_index
    }

    /// Returns the planned scratch size when checked arithmetic could
    /// represent it.
    pub fn planned_scratch_bytes(self) -> Option<u64> {
        self.planned_scratch_bytes
    }

    /// Returns the query scratch limit for a capacity failure.
    pub fn scratch_limit_bytes(self) -> Option<u64> {
        self.scratch_limit_bytes
    }

    fn with_point_index(mut self, point_index: usize) -> Self {
        self.point_index = Some(point_index);
        self
    }

    fn from_sample_failure(failure: QuerySampleFailure) -> Self {
        Self {
            reason: match failure {
                QuerySampleFailure::NonFiniteResult => QueryErrorReason::NonFiniteResult,
                QuerySampleFailure::NumericalIndeterminate => {
                    QueryErrorReason::NumericalIndeterminate
                }
            },
            point_index: None,
            planned_scratch_bytes: None,
            scratch_limit_bytes: None,
        }
    }

    fn from_planning_failure(failure: QueryPlanningFailure) -> Self {
        let (planned_scratch_bytes, scratch_limit_bytes) = match failure {
            QueryPlanningFailure::ArithmeticOverflow => (None, QUERY_SCRATCH_LIMIT_BYTES),
            QueryPlanningFailure::ScratchLimitExceeded {
                planned_scratch_bytes,
                limit_bytes,
            } => (Some(planned_scratch_bytes), limit_bytes),
        };
        Self {
            reason: QueryErrorReason::CapacityExceeded,
            point_index: None,
            planned_scratch_bytes,
            scratch_limit_bytes: Some(scratch_limit_bytes),
        }
    }
}

impl fmt::Display for QueryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.reason, self.point_index) {
            (QueryErrorReason::NonFiniteResult, Some(index)) => write!(
                formatter,
                "model evaluation at batch index {index} did not produce finite observables"
            ),
            (QueryErrorReason::NonFiniteResult, None) => {
                formatter.write_str("model evaluation did not produce finite observables")
            }
            (QueryErrorReason::NumericalIndeterminate, Some(index)) => write!(
                formatter,
                "model evaluation at batch index {index} remained numerically indeterminate"
            ),
            (QueryErrorReason::NumericalIndeterminate, None) => {
                formatter.write_str("model evaluation remained numerically indeterminate")
            }
            (QueryErrorReason::CapacityExceeded, _) => {
                formatter.write_str("logical query batch exceeds the query scratch capacity")
            }
        }
    }
}

impl Error for QueryError {}

#[cfg(test)]
mod tests {
    use super::{
        FieldSample, QUERY_SCRATCH_LIMIT_BYTES, QueryPlanningFailure, plan_query_capacity,
    };

    #[test]
    fn query_planning_admits_the_guarantee_and_rejects_the_first_over_limit_batch() {
        let guarantee = plan_query_capacity(100_000)
            .expect("100,000 value-gradient samples fit the query scratch envelope");
        assert!(guarantee.scratch_bytes <= QUERY_SCRATCH_LIMIT_BYTES);
        assert!(guarantee.chunk_len > 0);
        assert!(guarantee.chunk_scratch_bytes > 0);

        let bytes_per_sample = size_of::<FieldSample>() as u64;
        let bytes_per_recovered_sample = size_of::<super::RecoveredQueryResult>() as u64;
        assert_eq!(
            guarantee.scratch_bytes,
            100_000 * bytes_per_sample + guarantee.chunk_scratch_bytes
        );
        let maximum_points =
            ((QUERY_SCRATCH_LIMIT_BYTES - bytes_per_recovered_sample) / bytes_per_sample) as usize;
        let boundary = plan_query_capacity(maximum_points)
            .expect("the largest batch retaining one query-result slot remains admissible");
        assert!(boundary.scratch_bytes <= QUERY_SCRATCH_LIMIT_BYTES);

        let first_rejected_bytes =
            (maximum_points as u64 + 1) * bytes_per_sample + bytes_per_recovered_sample;
        assert!(first_rejected_bytes > QUERY_SCRATCH_LIMIT_BYTES);

        assert_eq!(
            plan_query_capacity(maximum_points + 1),
            Err(QueryPlanningFailure::ScratchLimitExceeded {
                planned_scratch_bytes: first_rejected_bytes,
                limit_bytes: QUERY_SCRATCH_LIMIT_BYTES,
            })
        );
    }

    #[test]
    fn query_capacity_failure_is_exposed_without_a_point_index() {
        let failure = plan_query_capacity(usize::MAX)
            .expect_err("unrepresentable batch scratch is rejected by checked planning");
        let error = super::QueryError::from_planning_failure(failure);

        assert_eq!(error.reason(), super::QueryErrorReason::CapacityExceeded);
        assert_eq!(error.point_index(), None);
        assert_eq!(error.planned_scratch_bytes(), None);
        assert_eq!(error.scratch_limit_bytes(), Some(QUERY_SCRATCH_LIMIT_BYTES));
    }
}
