//! Kernel configuration visible at the supported public fitting boundary.

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

/// Resolved dimensionless multiplier applied to the kernel's native FieldEnergy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FieldEnergyNormalization {
    factor: f64,
}

impl FieldEnergyNormalization {
    pub(crate) fn all_hard() -> Self {
        Self { factor: 1.0 }
    }

    /// Returns the finite positive resolved multiplier.
    pub fn factor(self) -> f64 {
        self.factor
    }
}
