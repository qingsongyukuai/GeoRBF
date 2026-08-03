//! Kernel configuration visible at the supported public fitting boundary.

use std::error::Error;
use std::fmt;

/// A resolved kernel kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum KernelKind {
    /// The scale-free three-dimensional Cubic kernel `r³` with complete Π₁.
    Cubic,
}

/// An immutable kernel configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KernelConfig {
    kind: KernelKind,
}

impl KernelConfig {
    /// Selects the supported Cubic Equality kernel contract.
    pub fn cubic() -> Self {
        Self {
            kind: KernelKind::Cubic,
        }
    }

    /// Returns the resolved kernel kind.
    pub fn kind(&self) -> KernelKind {
        self.kind
    }
}

impl Default for KernelConfig {
    fn default() -> Self {
        Self::cubic()
    }
}

/// Resolved physical multiplier applied to the kernel's native FieldEnergy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FieldEnergyNormalization {
    factor: f64,
}

impl FieldEnergyNormalization {
    /// Creates a finite, strictly positive FieldEnergy multiplier.
    pub fn try_new(factor: f64) -> Result<Self, FieldEnergyNormalizationError> {
        if !factor.is_finite() {
            return Err(FieldEnergyNormalizationError::NotFinite);
        }
        if factor <= 0.0 {
            return Err(FieldEnergyNormalizationError::NotPositive);
        }
        Ok(Self { factor })
    }

    pub(crate) fn all_hard() -> Self {
        Self { factor: 1.0 }
    }

    /// Returns the finite positive resolved multiplier.
    pub fn factor(self) -> f64 {
        self.factor
    }
}

/// A rejected FieldEnergy normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FieldEnergyNormalizationError {
    /// The normalization was NaN or infinite.
    NotFinite,
    /// The normalization was zero or negative.
    NotPositive,
}

impl fmt::Display for FieldEnergyNormalizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => formatter.write_str("FieldEnergy normalization is not finite"),
            Self::NotPositive => formatter.write_str("FieldEnergy normalization is not positive"),
        }
    }
}

impl Error for FieldEnergyNormalizationError {}
