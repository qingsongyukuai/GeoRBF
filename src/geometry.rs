//! Coordinates, frame metadata, unit labels, and global anisotropy.

use crate::cubic::{GlobalAnisotropyMetric as CubicMetric, MetricError};
use std::error::Error;
use std::fmt;

/// A finite point in the problem's declared input coordinate frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    components: [f64; 3],
}

impl Point3 {
    /// Creates a point, rejecting any non-finite coordinate.
    pub fn try_new(x: f64, y: f64, z: f64) -> Result<Self, GeometryError> {
        let components = [x, y, z];
        validate_components(components).map_err(|axis| GeometryError::NonFinitePoint { axis })?;
        Ok(Self { components })
    }

    /// Returns components in the declared frame's ordered basis.
    pub fn components(self) -> [f64; 3] {
        self.components
    }
}

/// A finite vector in the problem's declared input coordinate frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3 {
    components: [f64; 3],
}

impl Vector3 {
    /// Creates a vector, rejecting any non-finite component.
    pub fn try_new(x: f64, y: f64, z: f64) -> Result<Self, GeometryError> {
        let components = [x, y, z];
        validate_components(components).map_err(|axis| GeometryError::NonFiniteVector { axis })?;
        Ok(Self { components })
    }

    /// Returns components in the declared frame's ordered basis.
    pub fn components(self) -> [f64; 3] {
        self.components
    }
}

pub(crate) fn normalize_direction(direction: Vector3) -> Option<Vector3> {
    let components = direction.components();
    let scale = components
        .iter()
        .map(|component| component.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return None;
    }
    let scaled = components.map(|component| component / scale);
    let norm = scaled
        .into_iter()
        .map(|component| component * component)
        .sum::<f64>()
        .sqrt();
    let unit = scaled.map(|component| {
        let component = component / norm;
        if component == 0.0 { 0.0 } else { component }
    });
    Some(
        Vector3::try_new(unit[0], unit[1], unit[2])
            .expect("normalizing a finite nonzero vector produces finite components"),
    )
}

/// Handedness of the declared ordered orthogonal basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Handedness {
    /// A right-handed ordered basis.
    Right,
    /// A left-handed ordered basis.
    Left,
}

/// An opaque label for the common length unit of a problem.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LengthUnitLabel(Box<str>);

impl LengthUnitLabel {
    /// Owns the caller's unit label without interpreting or converting it.
    pub fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    /// Returns the caller-provided label.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for LengthUnitLabel {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// An opaque label for the scalar field unit of a problem.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FieldUnitLabel(Box<str>);

impl FieldUnitLabel {
    /// Owns the caller's unit label without interpreting or converting it.
    pub fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    /// Returns the caller-provided label.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for FieldUnitLabel {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Caller-declared semantics for the ordered coordinates used by a problem.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct InputCoordinateFrame {
    axis_labels: [Box<str>; 3],
    handedness: Handedness,
    length_unit: LengthUnitLabel,
}

impl InputCoordinateFrame {
    /// Creates a frame with three non-empty, distinct ordered axis labels.
    ///
    /// The labels declare an orthogonal Cartesian basis. GeoRBF does not infer
    /// axis meaning, a vertical direction, or a unit conversion from them.
    pub fn try_new<S>(
        axis_labels: [S; 3],
        handedness: Handedness,
        length_unit: LengthUnitLabel,
    ) -> Result<Self, GeometryError>
    where
        S: Into<Box<str>>,
    {
        let axis_labels = axis_labels.map(Into::into);
        for (axis, label) in axis_labels.iter().enumerate() {
            if label.is_empty() {
                return Err(GeometryError::EmptyAxisLabel { axis });
            }
        }
        for first in 0..3 {
            for second in (first + 1)..3 {
                if axis_labels[first] == axis_labels[second] {
                    return Err(GeometryError::DuplicateAxisLabel { first, second });
                }
            }
        }
        Ok(Self {
            axis_labels,
            handedness,
            length_unit,
        })
    }

    /// Returns axis labels in component order.
    pub fn axis_labels(&self) -> [&str; 3] {
        self.axis_labels.each_ref().map(|label| label.as_ref())
    }

    /// Returns the declared handedness.
    pub fn handedness(&self) -> Handedness {
        self.handedness
    }

    /// Returns the common length-unit label.
    pub fn length_unit(&self) -> &LengthUnitLabel {
        &self.length_unit
    }
}

/// A finite symmetric-positive-definite, determinant-one kernel metric.
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalAnisotropyMetric {
    inner: CubicMetric,
}

impl GlobalAnisotropyMetric {
    /// Returns the explicit identity metric used when no metric is configured.
    pub fn identity() -> Self {
        Self {
            inner: CubicMetric::identity(),
        }
    }

    /// Checks a matrix without symmetrizing, normalizing, or repairing it.
    pub fn try_from_matrix(matrix: [[f64; 3]; 3]) -> Result<Self, GeometryError> {
        CubicMetric::new(matrix)
            .map(|inner| Self { inner })
            .map_err(|reason| {
                GeometryError::InvalidGlobalAnisotropyMetric(match reason {
                    MetricError::NonFinite { row, column } => {
                        GlobalAnisotropyMetricError::NonFinite { row, column }
                    }
                    MetricError::NotSymmetric { row, column } => {
                        GlobalAnisotropyMetricError::NotSymmetric { row, column }
                    }
                    MetricError::NotPositiveDefinite => {
                        GlobalAnisotropyMetricError::NotPositiveDefinite
                    }
                    MetricError::NonFiniteDeterminant => {
                        GlobalAnisotropyMetricError::NonFiniteDeterminant
                    }
                    MetricError::DeterminantNotOne { determinant } => {
                        GlobalAnisotropyMetricError::DeterminantNotOne { determinant }
                    }
                })
            })
    }

    /// Returns the checked metric matrix.
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        self.inner.matrix()
    }

    pub(crate) fn as_cubic_metric(&self) -> CubicMetric {
        self.inner.clone()
    }
}

/// A rejected geometry or frame value.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum GeometryError {
    /// A point coordinate was not finite.
    NonFinitePoint { axis: usize },
    /// A vector component was not finite.
    NonFiniteVector { axis: usize },
    /// A frame axis label was empty.
    EmptyAxisLabel { axis: usize },
    /// Two frame axes had the same label.
    DuplicateAxisLabel { first: usize, second: usize },
    /// A global anisotropy metric violated its mathematical contract.
    InvalidGlobalAnisotropyMetric(GlobalAnisotropyMetricError),
}

/// Why a global anisotropy metric failed its checked constructor.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum GlobalAnisotropyMetricError {
    /// A matrix entry was not finite.
    NonFinite { row: usize, column: usize },
    /// Exact symmetry was violated.
    NotSymmetric { row: usize, column: usize },
    /// The matrix was not strictly positive definite.
    NotPositiveDefinite,
    /// Computing the determinant produced a non-finite value.
    NonFiniteDeterminant,
    /// The determinant was not one under the versioned policy.
    DeterminantNotOne { determinant: f64 },
}

impl fmt::Display for GeometryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinitePoint { axis } => {
                write!(formatter, "point coordinate {axis} is not finite")
            }
            Self::NonFiniteVector { axis } => {
                write!(formatter, "vector component {axis} is not finite")
            }
            Self::EmptyAxisLabel { axis } => write!(formatter, "frame axis {axis} has no label"),
            Self::DuplicateAxisLabel { first, second } => write!(
                formatter,
                "frame axes {first} and {second} have the same label"
            ),
            Self::InvalidGlobalAnisotropyMetric(reason) => {
                write!(formatter, "invalid global anisotropy metric: {reason:?}")
            }
        }
    }
}

impl Error for GeometryError {}

fn validate_components(components: [f64; 3]) -> Result<(), usize> {
    components
        .iter()
        .position(|component| !component.is_finite())
        .map_or(Ok(()), Err)
}
