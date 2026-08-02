//! Supported absolute field-value and complete-gradient observations.

use std::error::Error;
use std::fmt;

use crate::functional::SourceId;
use crate::geometry::{Point3, Vector3};
use crate::problem::{AddError, ProblemBuilder, ProblemInput, private};

/// A hard observation of an absolute scalar field value at one point.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldValueObservation {
    source_id: SourceId,
    location: Point3,
    value: f64,
}

impl FieldValueObservation {
    /// Creates a hard field-value observation.
    pub fn try_new(
        source_id: SourceId,
        location: Point3,
        value: f64,
    ) -> Result<Self, ObservationError> {
        if !value.is_finite() {
            return Err(ObservationError::NonFiniteFieldValue);
        }
        Ok(Self {
            source_id,
            location,
            value,
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

/// A rejected observation value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ObservationError {
    /// An absolute field value was NaN or infinite.
    NonFiniteFieldValue,
}

impl fmt::Display for ObservationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteFieldValue => formatter.write_str("field value is not finite"),
        }
    }
}

impl Error for ObservationError {}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ObservationInput {
    FieldValue(FieldValueObservation),
    Gradient(GradientObservation),
}

impl ObservationInput {
    pub(crate) fn source_id(&self) -> &SourceId {
        match self {
            Self::FieldValue(observation) => observation.source_id(),
            Self::Gradient(observation) => observation.source_id(),
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
