//! Supported field-value, complete-gradient, and tangent-direction observations.

use std::collections::BTreeSet;
use std::error::Error;
use std::fmt;

use crate::functional::{GroupId, SourceId};
use crate::geometry::{Point3, Vector3};
use crate::math::canonical_zero;
use crate::problem::{AddError, ProblemBuilder, ProblemInput, private};

/// A finite, exactly symmetric, strictly positive-definite covariance matrix.
///
/// The matrix is crate-owned and dynamically sized so the same checked value
/// can describe one vector observation or an explicitly ordered covariance
/// group. Entries are stored in row-major order and retain the caller's
/// physical residual units squared.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceMatrix {
    dimension: usize,
    entries: Box<[f64]>,
}

impl CovarianceMatrix {
    /// Checks a square covariance supplied as a fixed-size Rust array.
    pub fn try_new<const N: usize>(matrix: [[f64; N]; N]) -> Result<Self, CovarianceMatrixError> {
        Self::try_from_rows(matrix.into_iter().map(Vec::from).collect())
    }

    /// Checks a dynamically sized square covariance without repairing it.
    pub fn try_from_rows(rows: Vec<Vec<f64>>) -> Result<Self, CovarianceMatrixError> {
        let dimension = rows.len();
        if dimension == 0 {
            return Err(CovarianceMatrixError::Empty);
        }
        if let Some((row, actual)) = rows
            .iter()
            .enumerate()
            .find_map(|(row, values)| (values.len() != dimension).then_some((row, values.len())))
        {
            return Err(CovarianceMatrixError::NotSquare {
                row,
                expected: dimension,
                actual,
            });
        }
        for row in 0..dimension {
            for column in 0..dimension {
                if !rows[row][column].is_finite() {
                    return Err(CovarianceMatrixError::NonFinite { row, column });
                }
                if column > row && rows[row][column] != rows[column][row] {
                    return Err(CovarianceMatrixError::NotSymmetric { row, column });
                }
            }
        }
        let scale = rows
            .iter()
            .flatten()
            .map(|entry| entry.abs())
            .fold(0.0_f64, f64::max);
        if scale == 0.0 || !scaled_cholesky_is_positive_definite(&rows, scale) {
            return Err(CovarianceMatrixError::NotPositiveDefinite);
        }
        Ok(Self {
            dimension,
            entries: rows.into_iter().flatten().collect(),
        })
    }

    /// Returns the number of scalar residual components described.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns one row-major matrix entry, or `None` for an invalid index.
    pub fn entry(&self, row: usize, column: usize) -> Option<f64> {
        (row < self.dimension && column < self.dimension)
            .then(|| self.entries[row * self.dimension + column])
    }

    pub(crate) fn entries(&self) -> &[f64] {
        &self.entries
    }
}

fn scaled_cholesky_is_positive_definite(rows: &[Vec<f64>], scale: f64) -> bool {
    let dimension = rows.len();
    let mut lower = vec![0.0; dimension * dimension];
    for row in 0..dimension {
        for column in 0..=row {
            let product = (0..column)
                .map(|index| lower[row * dimension + index] * lower[column * dimension + index])
                .sum::<f64>();
            let remainder = rows[row][column] / scale - product;
            if row == column {
                if !remainder.is_finite() || remainder <= 0.0 {
                    return false;
                }
                lower[row * dimension + column] = remainder.sqrt();
            } else {
                let value = remainder / lower[column * dimension + column];
                if !value.is_finite() {
                    return false;
                }
                lower[row * dimension + column] = value;
            }
        }
    }
    true
}

/// Why a covariance matrix failed its checked constructor.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CovarianceMatrixError {
    /// A covariance must describe at least one residual component.
    Empty,
    /// A dynamic row did not match the matrix dimension.
    NotSquare {
        row: usize,
        expected: usize,
        actual: usize,
    },
    /// A matrix entry was NaN or infinite.
    NonFinite { row: usize, column: usize },
    /// Exact symmetry was violated.
    NotSymmetric { row: usize, column: usize },
    /// The matrix was not strictly positive definite.
    NotPositiveDefinite,
}

impl fmt::Display for CovarianceMatrixError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => formatter.write_str("a covariance matrix cannot be empty"),
            Self::NotSquare {
                row,
                expected,
                actual,
            } => write!(
                formatter,
                "covariance row {row} has length {actual}, expected {expected}"
            ),
            Self::NonFinite { row, column } => {
                write!(
                    formatter,
                    "covariance entry ({row}, {column}) is not finite"
                )
            }
            Self::NotSymmetric { row, column } => {
                write!(
                    formatter,
                    "covariance entries ({row}, {column}) and ({column}, {row}) differ"
                )
            }
            Self::NotPositiveDefinite => {
                formatter.write_str("covariance matrix is not strictly positive definite")
            }
        }
    }
}

impl Error for CovarianceMatrixError {}

/// The physical dimension shared by every member of a covariance group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CovarianceResidualDimension {
    /// Scalar field-value residuals.
    FieldValue,
    /// Gradient or directional-derivative residuals.
    FieldValuePerLength,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CovarianceGroupMember {
    FieldValue(FieldValueObservation),
    Gradient(GradientObservation),
    Tangent(TangentDirectionObservation),
}

impl CovarianceGroupMember {
    pub(crate) fn source_id(&self) -> &SourceId {
        match self {
            Self::FieldValue(observation) => observation.source_id(),
            Self::Gradient(observation) => observation.source_id(),
            Self::Tangent(observation) => observation.source_id(),
        }
    }

    pub(crate) fn dimension(&self) -> CovarianceResidualDimension {
        match self {
            Self::FieldValue(_) => CovarianceResidualDimension::FieldValue,
            Self::Gradient(_) | Self::Tangent(_) => {
                CovarianceResidualDimension::FieldValuePerLength
            }
        }
    }

    pub(crate) fn scalar_residual_count(&self) -> usize {
        match self {
            Self::Gradient(_) => 3,
            Self::FieldValue(_) | Self::Tangent(_) => 1,
        }
    }
}

/// Atomically constructs one named statistical covariance group.
#[derive(Debug)]
pub struct CovarianceGroupBuilder {
    group_id: GroupId,
    members: Vec<CovarianceGroupMember>,
    source_ids: BTreeSet<SourceId>,
    dimension: Option<CovarianceResidualDimension>,
    scalar_residual_count: usize,
}

impl CovarianceGroupBuilder {
    /// Starts an incomplete covariance group with a stable caller identity.
    pub fn new(group_id: GroupId) -> Self {
        Self {
            group_id,
            members: Vec::new(),
            source_ids: BTreeSet::new(),
            dimension: None,
            scalar_residual_count: 0,
        }
    }

    /// Adds one scalar field-value residual member.
    pub fn add_field_value_member(
        &mut self,
        source_id: SourceId,
        location: Point3,
        value: f64,
    ) -> Result<(), CovarianceGroupMemberAddError> {
        let observation = FieldValueObservation::try_new(source_id, location, value).map_err(
            |error| match error {
                ObservationError::NonFiniteFieldValue => {
                    CovarianceGroupMemberAddError::NonFiniteFieldValue
                }
                _ => unreachable!("a field-value constructor has one local failure"),
            },
        )?;
        self.add_member(CovarianceGroupMember::FieldValue(observation))
    }

    /// Adds one complete three-component Gradient residual member.
    pub fn add_gradient_member(
        &mut self,
        source_id: SourceId,
        location: Point3,
        gradient: Vector3,
    ) -> Result<(), CovarianceGroupMemberAddError> {
        self.add_member(CovarianceGroupMember::Gradient(GradientObservation::new(
            source_id, location, gradient,
        )))
    }

    /// Adds one scalar zero directional-derivative Tangent residual member.
    pub fn add_tangent_member(
        &mut self,
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
    ) -> Result<(), CovarianceGroupMemberAddError> {
        let observation = TangentDirectionObservation::try_new(source_id, location, direction)
            .map_err(|error| match error {
                ObservationError::ZeroTangentDirection => {
                    CovarianceGroupMemberAddError::ZeroTangentDirection
                }
                _ => unreachable!("a checked tangent member has one local failure"),
            })?;
        self.add_member(CovarianceGroupMember::Tangent(observation))
    }

    fn add_member(
        &mut self,
        member: CovarianceGroupMember,
    ) -> Result<(), CovarianceGroupMemberAddError> {
        let actual = member.dimension();
        if let Some(expected) = self.dimension {
            if actual != expected {
                return Err(CovarianceGroupMemberAddError::DimensionMismatch { expected, actual });
            }
        }
        let source_id = member.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(CovarianceGroupMemberAddError::DuplicateSourceId { source_id });
        }
        self.scalar_residual_count += member.scalar_residual_count();
        self.dimension = Some(actual);
        self.source_ids.insert(source_id);
        self.members.push(member);
        Ok(())
    }

    /// Finishes the non-empty group when covariance and residual dimensions match.
    pub fn build(
        self,
        covariance: CovarianceMatrix,
    ) -> Result<CovarianceGroup, CovarianceGroupBuildError> {
        if self.members.is_empty() {
            return Err(CovarianceGroupBuildError::EmptyGroup);
        }
        if covariance.dimension() != self.scalar_residual_count {
            return Err(CovarianceGroupBuildError::CovarianceDimensionMismatch {
                expected: self.scalar_residual_count,
                actual: covariance.dimension(),
            });
        }
        Ok(CovarianceGroup {
            group_id: self.group_id,
            members: self.members,
            dimension: self
                .dimension
                .expect("a non-empty covariance group has a physical dimension"),
            scalar_residual_count: self.scalar_residual_count,
            covariance,
        })
    }
}

/// A complete, immutable named covariance group.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceGroup {
    group_id: GroupId,
    members: Vec<CovarianceGroupMember>,
    dimension: CovarianceResidualDimension,
    scalar_residual_count: usize,
    covariance: CovarianceMatrix,
}

impl CovarianceGroup {
    /// Returns the stable caller-owned group identity.
    pub fn group_id(&self) -> &GroupId {
        &self.group_id
    }

    /// Returns the explicit member count, without flattening vector members.
    pub fn member_count(&self) -> usize {
        self.members.len()
    }

    /// Returns the flattened ordered residual dimension.
    pub fn scalar_residual_count(&self) -> usize {
        self.scalar_residual_count
    }

    /// Returns the common physical residual dimension.
    pub fn residual_dimension(&self) -> CovarianceResidualDimension {
        self.dimension
    }

    /// Returns the checked covariance in flattened component order.
    pub fn covariance(&self) -> &CovarianceMatrix {
        &self.covariance
    }

    pub(crate) fn members(&self) -> &[CovarianceGroupMember] {
        &self.members
    }
}

/// A rejected covariance-group draft mutation.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CovarianceGroupMemberAddError {
    /// A member reused a SourceId already owned by this draft.
    DuplicateSourceId { source_id: SourceId },
    /// A member's physical residual dimension differed from prior members.
    DimensionMismatch {
        expected: CovarianceResidualDimension,
        actual: CovarianceResidualDimension,
    },
    /// A field-value member target was NaN or infinite.
    NonFiniteFieldValue,
    /// A Tangent member supplied the zero direction.
    ZeroTangentDirection,
}

impl fmt::Display for CovarianceGroupMemberAddError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateSourceId { source_id } => {
                write!(
                    formatter,
                    "duplicate covariance member SourceId `{source_id}`"
                )
            }
            Self::DimensionMismatch { expected, actual } => write!(
                formatter,
                "covariance member dimension {actual:?} does not match {expected:?}"
            ),
            Self::NonFiniteFieldValue => formatter.write_str("field value is not finite"),
            Self::ZeroTangentDirection => formatter.write_str("tangent direction is zero"),
        }
    }
}

impl Error for CovarianceGroupMemberAddError {}

/// A rejected attempt to finish a covariance group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CovarianceGroupBuildError {
    /// A named covariance group must contain at least one member.
    EmptyGroup,
    /// The covariance dimension did not match the flattened member residuals.
    CovarianceDimensionMismatch { expected: usize, actual: usize },
}

impl fmt::Display for CovarianceGroupBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGroup => formatter.write_str("a covariance group cannot be empty"),
            Self::CovarianceDimensionMismatch { expected, actual } => write!(
                formatter,
                "covariance dimension {actual} does not match group residual dimension {expected}"
            ),
        }
    }
}

impl Error for CovarianceGroupBuildError {}

/// A non-statistical weight for one scalar or Euclidean-vector residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticPenalty {
    weight: f64,
}

impl QuadraticPenalty {
    /// Creates a finite, strictly positive quadratic penalty weight.
    pub fn try_new(weight: f64) -> Result<Self, QuadraticPenaltyError> {
        if !weight.is_finite() {
            return Err(QuadraticPenaltyError::NotFinite);
        }
        if weight <= 0.0 {
            return Err(QuadraticPenaltyError::NotPositive);
        }
        Ok(Self { weight })
    }

    /// Returns the penalty weight in inverse squared residual units.
    pub fn weight(self) -> f64 {
        self.weight
    }
}

/// A rejected quadratic-penalty weight.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum QuadraticPenaltyError {
    /// The weight was NaN or infinite.
    NotFinite,
    /// The weight was zero or negative.
    NotPositive,
}

impl fmt::Display for QuadraticPenaltyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => formatter.write_str("quadratic penalty is not finite"),
            Self::NotPositive => formatter.write_str("quadratic penalty is not positive"),
        }
    }
}

impl Error for QuadraticPenaltyError {}

/// A statistical standard deviation for one scalar or isotropic vector residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StandardDeviation {
    value: f64,
}

impl StandardDeviation {
    /// Creates a finite, strictly positive standard deviation.
    pub fn try_new(value: f64) -> Result<Self, StandardDeviationError> {
        if !value.is_finite() {
            return Err(StandardDeviationError::NotFinite);
        }
        if value <= 0.0 {
            return Err(StandardDeviationError::NotPositive);
        }
        Ok(Self { value })
    }

    /// Returns the standard deviation in the configured residual's physical units.
    pub fn value(self) -> f64 {
        self.value
    }
}

/// A rejected statistical standard deviation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum StandardDeviationError {
    /// The standard deviation was NaN or infinite.
    NotFinite,
    /// The standard deviation was zero or negative.
    NotPositive,
}

impl fmt::Display for StandardDeviationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => formatter.write_str("standard deviation is not finite"),
            Self::NotPositive => formatter.write_str("standard deviation is not positive"),
        }
    }
}

impl Error for StandardDeviationError {}

/// A hard or explicitly weighted soft observation of an absolute field value.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldValueObservation {
    source_id: SourceId,
    location: Point3,
    value: f64,
    configuration: FieldValueConfiguration,
}

impl FieldValueObservation {
    /// Creates a hard field-value observation.
    pub fn try_new(
        source_id: SourceId,
        location: Point3,
        value: f64,
    ) -> Result<Self, ObservationError> {
        Self::new(source_id, location, value, FieldValueConfiguration::Hard)
    }

    /// Creates a soft field-value observation with a non-statistical
    /// quadratic penalty.
    pub fn try_with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        value: f64,
        penalty: QuadraticPenalty,
    ) -> Result<Self, ObservationError> {
        Self::new(
            source_id,
            location,
            value,
            FieldValueConfiguration::QuadraticPenalty(penalty),
        )
    }

    /// Creates a soft field-value observation with statistical uncertainty.
    pub fn try_with_standard_deviation(
        source_id: SourceId,
        location: Point3,
        value: f64,
        standard_deviation: StandardDeviation,
    ) -> Result<Self, ObservationError> {
        Self::new(
            source_id,
            location,
            value,
            FieldValueConfiguration::StandardDeviation(standard_deviation),
        )
    }

    fn new(
        source_id: SourceId,
        location: Point3,
        value: f64,
        configuration: FieldValueConfiguration,
    ) -> Result<Self, ObservationError> {
        if !value.is_finite() {
            return Err(ObservationError::NonFiniteFieldValue);
        }
        Ok(Self {
            source_id,
            location,
            value,
            configuration,
        })
    }

    /// Returns the stable caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the observation location.
    pub fn location(&self) -> Point3 {
        self.location
    }

    /// Returns the observed absolute field value.
    pub fn value(&self) -> f64 {
        self.value
    }

    pub(crate) fn configuration(&self) -> FieldValueConfiguration {
        self.configuration
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum FieldValueConfiguration {
    Hard,
    QuadraticPenalty(QuadraticPenalty),
    StandardDeviation(StandardDeviation),
}

impl FieldValueConfiguration {
    pub(crate) fn is_soft(self) -> bool {
        !matches!(self, Self::Hard)
    }
}

/// A hard observation of the complete field gradient at one point.
#[derive(Debug, Clone, PartialEq)]
pub struct GradientObservation {
    source_id: SourceId,
    location: Point3,
    gradient: Vector3,
    configuration: GradientConfiguration,
}

impl GradientObservation {
    /// Creates a hard complete-gradient observation.
    ///
    /// `Point3` and `Vector3` have already checked finiteness, so this
    /// constructor cannot admit a partial or non-finite observation.
    pub fn new(source_id: SourceId, location: Point3, gradient: Vector3) -> Self {
        Self {
            source_id,
            location,
            gradient,
            configuration: GradientConfiguration::Hard,
        }
    }

    /// Creates a soft complete-gradient observation with one rotationally
    /// invariant Euclidean vector quadratic penalty.
    pub fn with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        gradient: Vector3,
        penalty: QuadraticPenalty,
    ) -> Self {
        Self {
            source_id,
            location,
            gradient,
            configuration: GradientConfiguration::QuadraticPenalty(penalty),
        }
    }

    /// Creates an isotropic statistical soft complete-gradient observation.
    pub fn with_standard_deviation(
        source_id: SourceId,
        location: Point3,
        gradient: Vector3,
        standard_deviation: StandardDeviation,
    ) -> Self {
        Self {
            source_id,
            location,
            gradient,
            configuration: GradientConfiguration::StandardDeviation(standard_deviation),
        }
    }

    /// Creates a statistical soft complete-gradient observation.
    pub fn try_with_covariance(
        source_id: SourceId,
        location: Point3,
        gradient: Vector3,
        covariance: CovarianceMatrix,
    ) -> Result<Self, ObservationError> {
        if covariance.dimension() != 3 {
            return Err(ObservationError::CovarianceDimensionMismatch {
                expected: 3,
                actual: covariance.dimension(),
            });
        }
        Ok(Self {
            source_id,
            location,
            gradient,
            configuration: GradientConfiguration::Covariance(covariance),
        })
    }

    /// Returns the stable caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the observation location.
    pub fn location(&self) -> Point3 {
        self.location
    }

    /// Returns the observed complete gradient in the declared input frame.
    pub fn gradient(&self) -> Vector3 {
        self.gradient
    }

    pub(crate) fn configuration(&self) -> &GradientConfiguration {
        &self.configuration
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum GradientConfiguration {
    Hard,
    QuadraticPenalty(QuadraticPenalty),
    StandardDeviation(StandardDeviation),
    Covariance(CovarianceMatrix),
}

impl GradientConfiguration {
    pub(crate) fn is_soft(&self) -> bool {
        !matches!(self, Self::Hard)
    }
}

/// A hard observation of an unoriented tangent direction at one point.
///
/// A tangent direction constrains only the directional derivative to zero. It
/// permits a zero gradient and therefore does not assert a regular level set,
/// a normal polarity, or a gradient magnitude. Use [`GradientObservation`]
/// when the complete gradient vector is observed. Normal-direction semantics
/// would additionally constrain gradient alignment, polarity or axial
/// equivalence, and nonzero slope; no normal observation is public in this
/// milestone.
#[derive(Debug, Clone, PartialEq)]
pub struct TangentDirectionObservation {
    source_id: SourceId,
    location: Point3,
    direction: Vector3,
    configuration: TangentConfiguration,
}

impl TangentDirectionObservation {
    /// Creates a hard, unoriented unit tangent-direction observation.
    ///
    /// The direction is normalized in the problem's physical input
    /// coordinates. Opposite and nonzero scaled vectors produce the same
    /// canonical axial representative. Non-finite vectors cannot be created
    /// by [`Vector3`], and a zero vector is rejected here.
    pub fn try_new(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
    ) -> Result<Self, ObservationError> {
        Self::new(source_id, location, direction, TangentConfiguration::Hard)
    }

    /// Creates a soft directional-derivative residual with a positive
    /// non-statistical quadratic penalty.
    pub fn try_with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        penalty: QuadraticPenalty,
    ) -> Result<Self, ObservationError> {
        Self::new(
            source_id,
            location,
            direction,
            TangentConfiguration::QuadraticPenalty(penalty),
        )
    }

    /// Creates a statistical soft directional-derivative residual.
    pub fn try_with_standard_deviation(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        standard_deviation: StandardDeviation,
    ) -> Result<Self, ObservationError> {
        Self::new(
            source_id,
            location,
            direction,
            TangentConfiguration::StandardDeviation(standard_deviation),
        )
    }

    fn new(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        configuration: TangentConfiguration,
    ) -> Result<Self, ObservationError> {
        let components = direction.components();
        let scale = components
            .iter()
            .map(|component| component.abs())
            .fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(ObservationError::ZeroTangentDirection);
        }
        let scaled = components.map(|component| component / scale);
        let norm = scaled
            .into_iter()
            .map(|component| component * component)
            .sum::<f64>()
            .sqrt();
        let mut unit = scaled.map(|component| canonical_zero(component / norm));
        if unit
            .iter()
            .find(|component| **component != 0.0)
            .is_some_and(|component| component.is_sign_negative())
        {
            unit = unit.map(|component| canonical_zero(-component));
        }
        let direction = Vector3::try_new(unit[0], unit[1], unit[2])
            .expect("normalizing a finite nonzero vector produces finite components");
        Ok(Self {
            source_id,
            location,
            direction,
            configuration,
        })
    }

    /// Returns the stable caller-owned source identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the observation location.
    pub fn location(&self) -> Point3 {
        self.location
    }

    /// Returns the canonical unit representative of the unoriented axis.
    pub fn direction(&self) -> Vector3 {
        self.direction
    }

    pub(crate) fn configuration(&self) -> TangentConfiguration {
        self.configuration
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum TangentConfiguration {
    Hard,
    QuadraticPenalty(QuadraticPenalty),
    StandardDeviation(StandardDeviation),
}

impl TangentConfiguration {
    pub(crate) fn is_soft(self) -> bool {
        !matches!(self, Self::Hard)
    }
}

/// A rejected observation value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ObservationError {
    /// An absolute field value was NaN or infinite.
    NonFiniteFieldValue,
    /// A tangent direction had no orientation because every component was zero.
    ZeroTangentDirection,
    /// A covariance did not match its residual block's scalar dimension.
    CovarianceDimensionMismatch { expected: usize, actual: usize },
}

impl fmt::Display for ObservationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteFieldValue => formatter.write_str("field value is not finite"),
            Self::ZeroTangentDirection => formatter.write_str("tangent direction is zero"),
            Self::CovarianceDimensionMismatch { expected, actual } => write!(
                formatter,
                "covariance dimension {actual} does not match residual dimension {expected}"
            ),
        }
    }
}

impl Error for ObservationError {}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ObservationInput {
    FieldValue(FieldValueObservation),
    Gradient(GradientObservation),
    TangentDirection(TangentDirectionObservation),
}

impl ObservationInput {
    pub(crate) fn source_id(&self) -> &SourceId {
        match self {
            Self::FieldValue(observation) => observation.source_id(),
            Self::Gradient(observation) => observation.source_id(),
            Self::TangentDirection(observation) => observation.source_id(),
        }
    }
}

impl private::Sealed for FieldValueObservation {}

impl ProblemInput for FieldValueObservation {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_observation(ObservationInput::FieldValue(self))
    }
}

impl private::Sealed for GradientObservation {}

impl ProblemInput for GradientObservation {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_observation(ObservationInput::Gradient(self))
    }
}

impl private::Sealed for TangentDirectionObservation {}

impl ProblemInput for TangentDirectionObservation {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_observation(ObservationInput::TangentDirection(self))
    }
}

impl private::Sealed for CovarianceGroup {}

impl ProblemInput for CovarianceGroup {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_covariance_group(self)
    }
}
