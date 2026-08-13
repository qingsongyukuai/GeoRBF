//! Safe geometric values shared by geological constraints.
//!
//! Source: `surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`.

use std::fmt;

/// Safe-construction failures for constraint geometry.
///
/// Frozen Surfe does not validate these inputs and can retain infinities or
/// produce NaNs. GeoRBF rejects them instead of reproducing invalid state.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ConstraintError {
    /// A coordinate, scalar, angle, or vector component is NaN or infinite.
    NonFiniteInput,
    /// A planar normal has no direction.
    ZeroNormal,
    /// A tangent vector has no direction.
    ZeroTangent,
    /// `acos(normal.z)` would be outside its real-valued domain.
    NormalZOutOfRange,
    /// Strike/dip conversion produced a zero or non-finite normal.
    DegenerateOrientation,
    /// Public Surfe orientation input accepts only polarity codes 0 and 1.
    InvalidPolarity,
}

impl fmt::Display for ConstraintError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFiniteInput => "constraint input must be finite",
            Self::ZeroNormal => "planar normal must be non-zero",
            Self::ZeroTangent => "tangent vector must be non-zero",
            Self::NormalZOutOfRange => "planar normal z component must be in [-1, 1]",
            Self::DegenerateOrientation => "strike and dip produce a degenerate normal",
            Self::InvalidPolarity => "polarity must be 0 (upright) or 1 (overturned)",
        })
    }
}

impl std::error::Error for ConstraintError {}

pub(crate) fn require_finite(values: &[f64]) -> Result<(), ConstraintError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(ConstraintError::NonFiniteInput)
    }
}

pub(crate) fn is_zero_vector(vector: [f64; 3]) -> bool {
    vector.into_iter().all(|component| component == 0.0)
}

/// A four-coordinate Surfe point with initialized scalar and vector outputs.
///
/// The `c` coordinate is retained because frozen kernel radius calculations
/// include it. Constraint API entry points use its default value, zero.
#[derive(Clone, Debug)]
pub struct Point {
    x: f64,
    y: f64,
    z: f64,
    c: f64,
    scalar_field: f64,
    field_normal: [f64; 3],
}

impl Point {
    /// Construct a point with Surfe's default fourth coordinate, `c = 0`.
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self, ConstraintError> {
        Self::with_c(x, y, z, 0.0)
    }

    /// Construct a point with an explicit fourth coordinate.
    pub fn with_c(x: f64, y: f64, z: f64, c: f64) -> Result<Self, ConstraintError> {
        require_finite(&[x, y, z, c])?;
        Ok(Self {
            x,
            y,
            z,
            c,
            scalar_field: 0.0,
            field_normal: [0.0; 3],
        })
    }

    pub const fn x(&self) -> f64 {
        self.x
    }

    pub const fn y(&self) -> f64 {
        self.y
    }

    pub const fn z(&self) -> f64 {
        self.z
    }

    pub const fn c(&self) -> f64 {
        self.c
    }

    pub const fn position(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    pub const fn scalar_field(&self) -> f64 {
        self.scalar_field
    }

    pub const fn vector_field(&self) -> [f64; 3] {
        self.field_normal
    }

    pub const fn nx_interp(&self) -> f64 {
        self.field_normal[0]
    }

    pub const fn ny_interp(&self) -> f64 {
        self.field_normal[1]
    }

    pub const fn nz_interp(&self) -> f64 {
        self.field_normal[2]
    }

    pub fn set_c(&mut self, c: f64) -> Result<(), ConstraintError> {
        require_finite(&[c])?;
        self.c = c;
        Ok(())
    }

    pub fn set_scalar_field(&mut self, scalar_field: f64) -> Result<(), ConstraintError> {
        require_finite(&[scalar_field])?;
        self.scalar_field = scalar_field;
        Ok(())
    }

    pub fn set_vector_field(&mut self, nx: f64, ny: f64, nz: f64) -> Result<(), ConstraintError> {
        require_finite(&[nx, ny, nz])?;
        self.field_normal = [nx, ny, nz];
        Ok(())
    }
}
