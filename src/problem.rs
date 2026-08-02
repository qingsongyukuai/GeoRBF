//! Problem construction, immutable snapshots, identity, and fit configuration.

use std::collections::BTreeSet;
use std::error::Error;
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::Arc;

use crate::fit::{self, FitFailure, FitSuccess};
pub use crate::functional::{GroupId, SourceId};
use crate::geometry::{FieldUnitLabel, GlobalAnisotropyMetric, InputCoordinateFrame};
use crate::kernel::{FieldEnergyNormalization, KernelConfig};
use crate::numerical::NumericalPolicyId;
use crate::observation::ObservationInput;

/// Resource request for a synchronous fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ThreadBudget {
    /// Let GeoRBF resolve the request under its numerical policy.
    Automatic,
    /// Request an exact, non-zero thread count.
    Exact(NonZeroUsize),
}

impl Default for ThreadBudget {
    fn default() -> Self {
        Self::Automatic
    }
}

/// Complete public fit configuration stored in a problem snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FitConfiguration {
    numerical_policy: NumericalPolicyId,
    thread_budget: ThreadBudget,
}

impl FitConfiguration {
    /// Returns the resolved, versioned numerical policy identity.
    pub fn numerical_policy(self) -> NumericalPolicyId {
        self.numerical_policy
    }

    /// Returns the requested thread budget.
    pub fn thread_budget(self) -> ThreadBudget {
        self.thread_budget
    }

    /// Replaces the resource request while retaining the complete policy.
    pub fn with_thread_budget(mut self, thread_budget: ThreadBudget) -> Self {
        self.thread_budget = thread_budget;
        self
    }
}

impl Default for FitConfiguration {
    fn default() -> Self {
        Self {
            numerical_policy: NumericalPolicyId::georbf_v1(),
            thread_budget: ThreadBudget::Automatic,
        }
    }
}

/// Sealed marker implemented by domain inputs accepted by [`ProblemBuilder`].
#[allow(private_bounds)]
pub trait ProblemInput: private::Sealed {
    #[doc(hidden)]
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError>;
}

pub(crate) mod private {
    pub trait Sealed {}
}

/// The only mutable stage of a GeoRBF problem.
#[derive(Debug)]
pub struct ProblemBuilder {
    input_coordinate_frame: InputCoordinateFrame,
    field_unit: FieldUnitLabel,
    observations: Vec<ObservationInput>,
    source_ids: BTreeSet<SourceId>,
    global_anisotropy_metric: Option<GlobalAnisotropyMetric>,
    fit_configuration: FitConfiguration,
}

impl ProblemBuilder {
    /// Starts a problem with explicit coordinate-frame and field-unit semantics.
    pub fn new(input_coordinate_frame: InputCoordinateFrame, field_unit: FieldUnitLabel) -> Self {
        Self {
            input_coordinate_frame,
            field_unit,
            observations: Vec::new(),
            source_ids: BTreeSet::new(),
            global_anisotropy_metric: None,
            fit_configuration: FitConfiguration::default(),
        }
    }

    /// Atomically adds one crate-owned problem input.
    pub fn add<T: ProblemInput>(&mut self, input: T) -> Result<(), AddError> {
        input.add_to(self)
    }

    /// Sets the one global kernel metric accepted by this problem.
    pub fn set_global_anisotropy_metric(
        &mut self,
        metric: GlobalAnisotropyMetric,
    ) -> Result<(), BuilderConfigurationError> {
        if self.global_anisotropy_metric.is_some() {
            return Err(BuilderConfigurationError::GlobalAnisotropyMetricAlreadySet);
        }
        self.global_anisotropy_metric = Some(metric);
        Ok(())
    }

    /// Replaces the complete fit configuration before snapshot creation.
    pub fn set_fit_configuration(&mut self, configuration: FitConfiguration) {
        self.fit_configuration = configuration;
    }

    /// Validates cross-record state and creates an owning immutable snapshot.
    pub fn build(mut self) -> Result<ProblemSnapshot, BuildFailure> {
        if self.observations.is_empty() {
            return Err(BuildFailure {
                builder: self,
                errors: vec![BuildError::NoObservations],
            });
        }
        self.observations
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        let data = ProblemData {
            input_coordinate_frame: self.input_coordinate_frame,
            field_unit: self.field_unit,
            observations: self.observations,
            resolved_kernel: KernelConfig::default(),
            field_energy_normalization: FieldEnergyNormalization::all_hard(),
            global_anisotropy_metric: self
                .global_anisotropy_metric
                .unwrap_or_else(GlobalAnisotropyMetric::identity),
            fit_configuration: self.fit_configuration,
        };
        Ok(ProblemSnapshot {
            inner: Arc::new(data),
        })
    }

    pub(crate) fn add_observation(
        &mut self,
        observation: ObservationInput,
    ) -> Result<(), AddError> {
        let source_id = observation.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.observations.push(observation);
        Ok(())
    }
}

/// An owning, immutable problem snapshot.
#[derive(Debug, Clone)]
pub struct ProblemSnapshot {
    pub(crate) inner: Arc<ProblemData>,
}

impl ProblemSnapshot {
    /// Fits this snapshot synchronously using its stored configuration.
    pub fn fit(&self) -> Result<FitSuccess, FitFailure> {
        fit::fit_snapshot(self)
    }

    /// Returns the caller-declared input coordinate frame.
    pub fn input_coordinate_frame(&self) -> &InputCoordinateFrame {
        &self.inner.input_coordinate_frame
    }

    /// Returns the caller-declared scalar field unit label.
    pub fn field_unit(&self) -> &FieldUnitLabel {
        &self.inner.field_unit
    }

    /// Returns the resolved kernel configuration.
    pub fn resolved_kernel(&self) -> &KernelConfig {
        &self.inner.resolved_kernel
    }

    /// Returns the resolved all-hard FieldEnergy normalization.
    pub fn field_energy_normalization(&self) -> FieldEnergyNormalization {
        self.inner.field_energy_normalization
    }

    /// Returns the resolved global anisotropy metric, including identity.
    pub fn global_anisotropy_metric(&self) -> &GlobalAnisotropyMetric {
        &self.inner.global_anisotropy_metric
    }

    /// Returns the stored fit configuration.
    pub fn fit_configuration(&self) -> FitConfiguration {
        self.inner.fit_configuration
    }

    /// Returns the number of top-level observations.
    pub fn observation_count(&self) -> usize {
        self.inner.observations.len()
    }
}

#[derive(Debug)]
pub(crate) struct ProblemData {
    pub(crate) input_coordinate_frame: InputCoordinateFrame,
    pub(crate) field_unit: FieldUnitLabel,
    pub(crate) observations: Vec<ObservationInput>,
    pub(crate) resolved_kernel: KernelConfig,
    pub(crate) field_energy_normalization: FieldEnergyNormalization,
    pub(crate) global_anisotropy_metric: GlobalAnisotropyMetric,
    pub(crate) fit_configuration: FitConfiguration,
}

/// An atomically rejected builder insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum AddError {
    /// A top-level input reused an existing stable source identity.
    DuplicateSourceId { source_id: SourceId },
}

impl fmt::Display for AddError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateSourceId { source_id } => {
                write!(formatter, "duplicate SourceId `{source_id}`")
            }
        }
    }
}

impl Error for AddError {}

/// An atomically rejected problem-level configuration change.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BuilderConfigurationError {
    /// The problem already has its one global anisotropy metric.
    GlobalAnisotropyMetricAlreadySet,
}

impl fmt::Display for BuilderConfigurationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GlobalAnisotropyMetricAlreadySet => {
                formatter.write_str("a global anisotropy metric is already set")
            }
        }
    }
}

impl Error for BuilderConfigurationError {}

/// A deterministic cross-record build error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BuildError {
    /// The problem contains no observations.
    NoObservations,
}

/// Failed snapshot construction, retaining the original builder for repair.
#[derive(Debug)]
pub struct BuildFailure {
    builder: ProblemBuilder,
    errors: Vec<BuildError>,
}

impl BuildFailure {
    /// Returns all build errors in deterministic order.
    pub fn errors(&self) -> &[BuildError] {
        &self.errors
    }

    /// Recovers the unchanged builder for repair and retry.
    pub fn into_builder(self) -> ProblemBuilder {
        self.builder
    }
}

impl fmt::Display for BuildFailure {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "problem build failed with {} error(s)",
            self.errors.len()
        )
    }
}

impl Error for BuildFailure {}
