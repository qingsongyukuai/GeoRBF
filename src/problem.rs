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
use crate::observation::{CovarianceGroup, ObservationInput};
pub use crate::relation::StratigraphicFieldDirection;
use crate::relation::{
    AdditiveFieldGauge, AdditiveFieldGaugeReference, DirectionalDerivativeInterval,
    FieldSeparationInterval, FieldValueBound, PointToLevelSetRelation, SharedLevelSetInput,
    SharedLevelSetRelationInput,
};

/// Resource request for a synchronous fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ThreadBudget {
    /// Let GeoRBF resolve the request under its numerical policy.
    Automatic,
    /// Request an exact, non-zero thread count.
    ///
    /// The current Cubic Equality path admits `Exact(1)`; a different exact
    /// request is rejected by `ProblemBuilder::build` before fitting.
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
pub trait ProblemInput: private::Sealed + Sized {
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
    covariance_groups: Vec<CovarianceGroup>,
    shared_level_sets: Vec<SharedLevelSetInput>,
    additive_field_gauges: Vec<AdditiveFieldGauge>,
    field_value_bounds: Vec<FieldValueBound>,
    directional_derivative_intervals: Vec<DirectionalDerivativeInterval>,
    field_separation_intervals: Vec<FieldSeparationInterval>,
    point_to_level_set_relations: Vec<PointToLevelSetRelation>,
    shared_level_set_relations: Vec<SharedLevelSetRelationInput>,
    source_ids: BTreeSet<SourceId>,
    group_ids: BTreeSet<GroupId>,
    shared_level_group_ids: BTreeSet<GroupId>,
    field_energy_normalization: Option<FieldEnergyNormalization>,
    global_anisotropy_metric: Option<GlobalAnisotropyMetric>,
    stratigraphic_field_direction: Option<StratigraphicFieldDirection>,
    fit_configuration: FitConfiguration,
}

impl ProblemBuilder {
    /// Starts a problem with explicit coordinate-frame and field-unit semantics.
    pub fn new(input_coordinate_frame: InputCoordinateFrame, field_unit: FieldUnitLabel) -> Self {
        Self {
            input_coordinate_frame,
            field_unit,
            observations: Vec::new(),
            covariance_groups: Vec::new(),
            shared_level_sets: Vec::new(),
            additive_field_gauges: Vec::new(),
            field_value_bounds: Vec::new(),
            directional_derivative_intervals: Vec::new(),
            field_separation_intervals: Vec::new(),
            point_to_level_set_relations: Vec::new(),
            shared_level_set_relations: Vec::new(),
            source_ids: BTreeSet::new(),
            group_ids: BTreeSet::new(),
            shared_level_group_ids: BTreeSet::new(),
            field_energy_normalization: None,
            global_anisotropy_metric: None,
            stratigraphic_field_direction: None,
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

    /// Sets the physical scale between FieldEnergy and soft losses.
    pub fn set_field_energy_normalization(
        &mut self,
        normalization: FieldEnergyNormalization,
    ) -> Result<(), BuilderConfigurationError> {
        if self.field_energy_normalization.is_some() {
            return Err(BuilderConfigurationError::FieldEnergyNormalizationAlreadySet);
        }
        self.field_energy_normalization = Some(normalization);
        Ok(())
    }

    /// Sets the one explicit mapping from stratigraphic age to field value.
    pub fn set_stratigraphic_field_direction(
        &mut self,
        direction: StratigraphicFieldDirection,
    ) -> Result<(), BuilderConfigurationError> {
        if self.stratigraphic_field_direction.is_some() {
            return Err(BuilderConfigurationError::StratigraphicFieldDirectionAlreadySet);
        }
        self.stratigraphic_field_direction = Some(direction);
        Ok(())
    }

    /// Replaces the complete fit configuration before snapshot creation.
    pub fn set_fit_configuration(&mut self, configuration: FitConfiguration) {
        self.fit_configuration = configuration;
    }

    /// Validates cross-record state and creates an owning immutable snapshot.
    pub fn build(mut self) -> Result<ProblemSnapshot, BuildFailure> {
        let mut errors = Vec::new();
        if self.source_ids.is_empty() {
            errors.push(BuildError::NoObservations);
        }
        let has_soft_relation = !self.covariance_groups.is_empty()
            || self
                .shared_level_set_relations
                .iter()
                .any(SharedLevelSetRelationInput::is_soft)
            || self.field_value_bounds.iter().any(FieldValueBound::is_soft)
            || self
                .directional_derivative_intervals
                .iter()
                .any(DirectionalDerivativeInterval::is_soft)
            || self
                .field_separation_intervals
                .iter()
                .any(FieldSeparationInterval::is_soft)
            || self
                .point_to_level_set_relations
                .iter()
                .any(PointToLevelSetRelation::is_soft)
            || self
                .observations
                .iter()
                .any(|observation| match observation {
                    ObservationInput::FieldValue(value) => value.configuration().is_soft(),
                    ObservationInput::Gradient(gradient) => gradient.configuration().is_soft(),
                    ObservationInput::TangentDirection(tangent) => {
                        tangent.configuration().is_soft()
                    }
                });
        if has_soft_relation && self.field_energy_normalization.is_none() {
            errors.push(BuildError::MissingFieldEnergyNormalization);
        }
        if self.stratigraphic_field_direction.is_none()
            && self
                .shared_level_set_relations
                .iter()
                .any(SharedLevelSetRelationInput::is_stratigraphic_age_relation)
        {
            errors.push(BuildError::MissingStratigraphicFieldDirection);
        }
        let mut dangling_references = self
            .additive_field_gauges
            .iter()
            .filter_map(|gauge| match gauge.reference() {
                AdditiveFieldGaugeReference::Point(_) => None,
                AdditiveFieldGaugeReference::LevelSet(group_id)
                    if !self.shared_level_group_ids.contains(group_id) =>
                {
                    Some((gauge.source_id().clone(), group_id.clone()))
                }
                AdditiveFieldGaugeReference::LevelSet(_) => None,
            })
            .chain(self.shared_level_set_relations.iter().flat_map(|relation| {
                relation
                    .declared_group_ids()
                    .into_iter()
                    .filter(|group_id| !self.shared_level_group_ids.contains(*group_id))
                    .map(|group_id| (relation.source_id().clone(), group_id.clone()))
            }))
            .chain(self.field_separation_intervals.iter().flat_map(|relation| {
                [relation.reference_group_id(), relation.target_group_id()]
                    .into_iter()
                    .filter(|group_id| !self.shared_level_group_ids.contains(*group_id))
                    .map(|group_id| (relation.source_id().clone(), group_id.clone()))
            }))
            .chain(
                self.point_to_level_set_relations
                    .iter()
                    .filter(|relation| !self.shared_level_group_ids.contains(relation.group_id()))
                    .map(|relation| (relation.source_id().clone(), relation.group_id().clone())),
            )
            .collect::<Vec<_>>();
        dangling_references.sort();
        dangling_references.dedup();
        errors.extend(
            dangling_references
                .into_iter()
                .map(|(source_id, group_id)| BuildError::UnknownGroupReference {
                    source_id,
                    group_id,
                }),
        );
        if let ThreadBudget::Exact(count) = self.fit_configuration.thread_budget {
            if count.get() != 1 {
                errors.push(BuildError::UnsupportedThreadBudget {
                    requested: count.get(),
                });
            }
        }
        if !errors.is_empty() {
            return Err(BuildFailure {
                builder: self,
                errors,
            });
        }
        self.observations
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        self.shared_level_sets
            .sort_by(|left, right| left.group_id().cmp(right.group_id()));
        self.covariance_groups
            .sort_by(|left, right| left.group_id().cmp(right.group_id()));
        self.additive_field_gauges
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        self.field_value_bounds
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        self.directional_derivative_intervals
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        self.field_separation_intervals
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        self.point_to_level_set_relations
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        self.shared_level_set_relations
            .sort_by(|left, right| left.source_id().cmp(right.source_id()));
        let data = ProblemData {
            input_coordinate_frame: self.input_coordinate_frame,
            field_unit: self.field_unit,
            observations: self.observations,
            covariance_groups: self.covariance_groups,
            shared_level_sets: self.shared_level_sets,
            additive_field_gauges: self.additive_field_gauges,
            field_value_bounds: self.field_value_bounds,
            directional_derivative_intervals: self.directional_derivative_intervals,
            field_separation_intervals: self.field_separation_intervals,
            point_to_level_set_relations: self.point_to_level_set_relations,
            shared_level_set_relations: self.shared_level_set_relations,
            source_count: self.source_ids.len(),
            resolved_kernel: KernelConfig::default(),
            field_energy_normalization: self
                .field_energy_normalization
                .unwrap_or_else(FieldEnergyNormalization::all_hard),
            global_anisotropy_metric: self
                .global_anisotropy_metric
                .unwrap_or_else(GlobalAnisotropyMetric::identity),
            stratigraphic_field_direction: self.stratigraphic_field_direction,
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

    pub(crate) fn add_shared_level_set(
        &mut self,
        group: SharedLevelSetInput,
    ) -> Result<(), AddError> {
        let group_id = group.group_id().clone();
        if self.group_ids.contains(&group_id) {
            return Err(AddError::DuplicateGroupId { group_id });
        }
        if let Some(member) = group
            .members()
            .iter()
            .find(|member| self.source_ids.contains(member.source_id()))
        {
            return Err(AddError::DuplicateSourceId {
                source_id: member.source_id().clone(),
            });
        }
        self.group_ids.insert(group_id);
        self.shared_level_group_ids.insert(group.group_id().clone());
        self.source_ids.extend(
            group
                .members()
                .iter()
                .map(|member| member.source_id().clone()),
        );
        self.shared_level_sets.push(group);
        Ok(())
    }

    pub(crate) fn add_covariance_group(&mut self, group: CovarianceGroup) -> Result<(), AddError> {
        let group_id = group.group_id().clone();
        if self.group_ids.contains(&group_id) {
            return Err(AddError::DuplicateGroupId { group_id });
        }
        if let Some(member) = group
            .members()
            .iter()
            .find(|member| self.source_ids.contains(member.source_id()))
        {
            return Err(AddError::DuplicateSourceId {
                source_id: member.source_id().clone(),
            });
        }
        self.group_ids.insert(group_id);
        self.source_ids.extend(
            group
                .members()
                .iter()
                .map(|member| member.source_id().clone()),
        );
        self.covariance_groups.push(group);
        Ok(())
    }

    pub(crate) fn add_additive_field_gauge(
        &mut self,
        gauge: AdditiveFieldGauge,
    ) -> Result<(), AddError> {
        let source_id = gauge.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.additive_field_gauges.push(gauge);
        Ok(())
    }

    pub(crate) fn add_field_value_bound(&mut self, bound: FieldValueBound) -> Result<(), AddError> {
        let source_id = bound.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.field_value_bounds.push(bound);
        Ok(())
    }

    pub(crate) fn add_directional_derivative_interval(
        &mut self,
        interval: DirectionalDerivativeInterval,
    ) -> Result<(), AddError> {
        let source_id = interval.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.directional_derivative_intervals.push(interval);
        Ok(())
    }

    pub(crate) fn add_shared_level_set_relation(
        &mut self,
        relation: SharedLevelSetRelationInput,
    ) -> Result<(), AddError> {
        let source_id = relation.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.shared_level_set_relations.push(relation);
        Ok(())
    }

    pub(crate) fn add_field_separation_interval(
        &mut self,
        interval: FieldSeparationInterval,
    ) -> Result<(), AddError> {
        let source_id = interval.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.field_separation_intervals.push(interval);
        Ok(())
    }

    pub(crate) fn add_point_to_level_set_relation(
        &mut self,
        relation: PointToLevelSetRelation,
    ) -> Result<(), AddError> {
        let source_id = relation.source_id().clone();
        if self.source_ids.contains(&source_id) {
            return Err(AddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id);
        self.point_to_level_set_relations.push(relation);
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

    /// Returns the resolved FieldEnergy normalization, including the all-hard default.
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

    /// Returns the number of general mathematical shared level sets.
    pub fn shared_level_set_count(&self) -> usize {
        self.inner
            .shared_level_sets
            .iter()
            .filter(|group| !group.is_horizon())
            .count()
    }

    /// Returns the number of geological horizons.
    pub fn horizon_count(&self) -> usize {
        self.inner
            .shared_level_sets
            .iter()
            .filter(|group| group.is_horizon())
            .count()
    }

    /// Returns the number of complete named covariance groups.
    pub fn covariance_group_count(&self) -> usize {
        self.inner.covariance_groups.len()
    }

    /// Returns the number of caller-owned Field Value Bound relations.
    pub fn field_value_bound_count(&self) -> usize {
        self.inner.field_value_bounds.len()
    }

    /// Returns the caller-owned Directional Derivative Interval count.
    pub fn directional_derivative_interval_count(&self) -> usize {
        self.inner.directional_derivative_intervals.len()
    }

    /// Returns the caller-owned Field Separation Interval count.
    pub fn field_separation_interval_count(&self) -> usize {
        self.inner.field_separation_intervals.len()
    }

    /// Returns the caller-owned Point to Level Set Relation count.
    pub fn point_to_level_set_relation_count(&self) -> usize {
        self.inner.point_to_level_set_relations.len()
    }

    /// Returns the explicitly configured stratigraphic field direction.
    pub fn stratigraphic_field_direction(&self) -> Option<StratigraphicFieldDirection> {
        self.inner.stratigraphic_field_direction
    }

    /// Returns the number of caller-owned Younger Than and Older Than relations.
    pub fn stratigraphic_age_relation_count(&self) -> usize {
        self.inner
            .shared_level_set_relations
            .iter()
            .filter(|relation| relation.is_stratigraphic_age_relation())
            .count()
    }

    /// Returns the number of caller-owned non-strict Field Level Order relations.
    pub fn field_level_order_count(&self) -> usize {
        self.inner
            .shared_level_set_relations
            .iter()
            .filter(|relation| relation.is_field_level_order())
            .count()
    }

    /// Returns the number of independently identified caller sources.
    pub fn source_count(&self) -> usize {
        self.inner.source_count
    }
}

#[derive(Debug)]
pub(crate) struct ProblemData {
    pub(crate) input_coordinate_frame: InputCoordinateFrame,
    pub(crate) field_unit: FieldUnitLabel,
    pub(crate) observations: Vec<ObservationInput>,
    pub(crate) covariance_groups: Vec<CovarianceGroup>,
    pub(crate) shared_level_sets: Vec<SharedLevelSetInput>,
    pub(crate) additive_field_gauges: Vec<AdditiveFieldGauge>,
    pub(crate) field_value_bounds: Vec<FieldValueBound>,
    pub(crate) directional_derivative_intervals: Vec<DirectionalDerivativeInterval>,
    pub(crate) field_separation_intervals: Vec<FieldSeparationInterval>,
    pub(crate) point_to_level_set_relations: Vec<PointToLevelSetRelation>,
    pub(crate) shared_level_set_relations: Vec<SharedLevelSetRelationInput>,
    pub(crate) source_count: usize,
    pub(crate) resolved_kernel: KernelConfig,
    pub(crate) field_energy_normalization: FieldEnergyNormalization,
    pub(crate) global_anisotropy_metric: GlobalAnisotropyMetric,
    pub(crate) stratigraphic_field_direction: Option<StratigraphicFieldDirection>,
    pub(crate) fit_configuration: FitConfiguration,
}

/// An atomically rejected builder insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum AddError {
    /// A top-level input reused an existing stable source identity.
    DuplicateSourceId { source_id: SourceId },
    /// A complete group reused an existing stable group identity.
    DuplicateGroupId { group_id: GroupId },
}

impl fmt::Display for AddError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateSourceId { source_id } => {
                write!(formatter, "duplicate SourceId `{source_id}`")
            }
            Self::DuplicateGroupId { group_id } => {
                write!(formatter, "duplicate GroupId `{group_id}`")
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
    /// The problem already has its one FieldEnergy normalization.
    FieldEnergyNormalizationAlreadySet,
    /// The one problem-level stratigraphic field direction is already configured.
    StratigraphicFieldDirectionAlreadySet,
}

impl fmt::Display for BuilderConfigurationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GlobalAnisotropyMetricAlreadySet => {
                formatter.write_str("a global anisotropy metric is already set")
            }
            Self::FieldEnergyNormalizationAlreadySet => {
                formatter.write_str("a FieldEnergy normalization is already set")
            }
            Self::StratigraphicFieldDirectionAlreadySet => {
                formatter.write_str("a stratigraphic field direction is already set")
            }
        }
    }
}

impl Error for BuilderConfigurationError {}

/// A deterministic cross-record build error.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum BuildError {
    /// The problem contains no observations.
    NoObservations,
    /// At least one soft relation exists without an explicit FieldEnergy scale.
    MissingFieldEnergyNormalization,
    /// An age relation exists without an explicit stratigraphic field direction.
    MissingStratigraphicFieldDirection,
    /// The current public Cubic Equality path is intentionally sequential.
    UnsupportedThreadBudget { requested: usize },
    /// A relation references a GroupId absent from the completed snapshot.
    UnknownGroupReference {
        source_id: SourceId,
        group_id: GroupId,
    },
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
