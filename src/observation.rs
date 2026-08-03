//! Supported field-value, complete-gradient, and tangent-direction observations.

use std::error::Error;
use std::fmt;

use crate::functional::SourceId;
use crate::geometry::{Point3, Vector3};
use crate::math::canonical_zero;
use crate::problem::{AddError, ProblemBuilder, ProblemInput, private};

/// A non-statistical weight for one quadratic field-value residual.
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

    /// Returns the penalty weight in inverse squared field-value units.
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

/// A statistical standard deviation for one scalar field-value residual.
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

    /// Returns the standard deviation in field-value units.
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
        }
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
}

/// A rejected observation value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ObservationError {
    /// An absolute field value was NaN or infinite.
    NonFiniteFieldValue,
    /// A tangent direction had no orientation because every component was zero.
    ZeroTangentDirection,
}

impl fmt::Display for ObservationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteFieldValue => formatter.write_str("field value is not finite"),
            Self::ZeroTangentDirection => formatter.write_str("tangent direction is zero"),
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
