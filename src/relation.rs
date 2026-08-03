//! Shared levels, gauges, and checked scalar affine field relations.

use std::collections::BTreeSet;
use std::error::Error;
use std::fmt;

use crate::functional::{GroupId, SemanticRolePath, SourceId};
use crate::geometry::{Point3, Vector3, normalize_direction};
use crate::observation::QuadraticPenalty;
use crate::problem::{AddError, ProblemBuilder, ProblemInput, private};

/// A finite positive weight for one nonnegative scalar violation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearViolationPenalty {
    weight: f64,
}

impl LinearViolationPenalty {
    /// Creates a finite, strictly positive linear violation penalty.
    pub fn try_new(weight: f64) -> Result<Self, LinearViolationPenaltyError> {
        if !weight.is_finite() {
            return Err(LinearViolationPenaltyError::NotFinite);
        }
        if weight <= 0.0 {
            return Err(LinearViolationPenaltyError::NotPositive);
        }
        Ok(Self { weight })
    }

    /// Returns the weight in inverse violation units.
    pub fn weight(self) -> f64 {
        self.weight
    }
}

/// A rejected linear violation penalty.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum LinearViolationPenaltyError {
    /// The weight was NaN or infinite.
    NotFinite,
    /// The weight was zero or negative.
    NotPositive,
}

impl fmt::Display for LinearViolationPenaltyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => formatter.write_str("linear violation penalty is not finite"),
            Self::NotPositive => formatter.write_str("linear violation penalty is not positive"),
        }
    }
}

impl Error for LinearViolationPenaltyError {}

/// One legal positive loss applied to a nonnegative field-value violation.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FieldValueViolationPenalty {
    /// A half-weighted squared loss, `1/2 * weight * violation^2`.
    Quadratic(QuadraticPenalty),
    /// A weighted absolute one-sided loss, `weight * violation`.
    Linear(LinearViolationPenalty),
}

/// One legal positive loss applied to a nonnegative derivative violation.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum DirectionalDerivativeViolationPenalty {
    /// A half-weighted squared loss, `1/2 * weight * violation^2`.
    Quadratic(QuadraticPenalty),
    /// A weighted absolute one-sided loss, `weight * violation`.
    Linear(LinearViolationPenalty),
}

/// One legal positive loss applied to a Field Separation Interval side.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FieldSeparationViolationPenalty {
    /// A half-weighted squared loss, `1/2 * weight * violation^2`.
    Quadratic(QuadraticPenalty),
    /// A weighted one-sided loss, `weight * violation`.
    Linear(LinearViolationPenalty),
}

/// A finite signed field-value interval between two ordered shared level sets.
///
/// The constrained quantity is `target - reference`. It has field-value units
/// and carries no age, spatial-side, distance, or physical-thickness meaning.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldSeparationInterval {
    source_id: SourceId,
    reference_group_id: GroupId,
    target_group_id: GroupId,
    lower: AffineBoundSide,
    upper: AffineBoundSide,
}

impl FieldSeparationInterval {
    /// Creates one hard closed interval for `target - reference`.
    pub fn try_hard(
        source_id: SourceId,
        reference_group_id: GroupId,
        target_group_id: GroupId,
        lower: f64,
        upper: f64,
    ) -> Result<Self, FieldSeparationIntervalError> {
        Self::new(
            source_id,
            reference_group_id,
            target_group_id,
            (lower, AffineBoundConfiguration::Hard),
            (upper, AffineBoundConfiguration::Hard),
        )
    }

    /// Creates a soft interval with independent quadratic penalties per side.
    #[allow(clippy::too_many_arguments)]
    pub fn try_with_quadratic_penalties(
        source_id: SourceId,
        reference_group_id: GroupId,
        target_group_id: GroupId,
        lower: f64,
        lower_penalty: QuadraticPenalty,
        upper: f64,
        upper_penalty: QuadraticPenalty,
    ) -> Result<Self, FieldSeparationIntervalError> {
        Self::new(
            source_id,
            reference_group_id,
            target_group_id,
            (
                lower,
                AffineBoundConfiguration::QuadraticPenalty(lower_penalty),
            ),
            (
                upper,
                AffineBoundConfiguration::QuadraticPenalty(upper_penalty),
            ),
        )
    }

    /// Creates a soft interval with independent linear penalties per side.
    #[allow(clippy::too_many_arguments)]
    pub fn try_with_linear_violation_penalties(
        source_id: SourceId,
        reference_group_id: GroupId,
        target_group_id: GroupId,
        lower: f64,
        lower_penalty: LinearViolationPenalty,
        upper: f64,
        upper_penalty: LinearViolationPenalty,
    ) -> Result<Self, FieldSeparationIntervalError> {
        Self::new(
            source_id,
            reference_group_id,
            target_group_id,
            (
                lower,
                AffineBoundConfiguration::LinearViolationPenalty(lower_penalty),
            ),
            (
                upper,
                AffineBoundConfiguration::LinearViolationPenalty(upper_penalty),
            ),
        )
    }

    /// Creates a soft interval whose sides independently select a legal loss.
    #[allow(clippy::too_many_arguments)]
    pub fn try_with_violation_penalties(
        source_id: SourceId,
        reference_group_id: GroupId,
        target_group_id: GroupId,
        lower: f64,
        lower_penalty: FieldSeparationViolationPenalty,
        upper: f64,
        upper_penalty: FieldSeparationViolationPenalty,
    ) -> Result<Self, FieldSeparationIntervalError> {
        let configuration = |penalty| match penalty {
            FieldSeparationViolationPenalty::Quadratic(penalty) => {
                AffineBoundConfiguration::QuadraticPenalty(penalty)
            }
            FieldSeparationViolationPenalty::Linear(penalty) => {
                AffineBoundConfiguration::LinearViolationPenalty(penalty)
            }
        };
        Self::new(
            source_id,
            reference_group_id,
            target_group_id,
            (lower, configuration(lower_penalty)),
            (upper, configuration(upper_penalty)),
        )
    }

    fn new(
        source_id: SourceId,
        reference_group_id: GroupId,
        target_group_id: GroupId,
        lower: (f64, AffineBoundConfiguration),
        upper: (f64, AffineBoundConfiguration),
    ) -> Result<Self, FieldSeparationIntervalError> {
        if reference_group_id == target_group_id {
            return Err(FieldSeparationIntervalError::SelfReference {
                group_id: reference_group_id,
            });
        }
        if !lower.0.is_finite() || !upper.0.is_finite() {
            return Err(FieldSeparationIntervalError::NonFiniteBound);
        }
        if lower.0 > upper.0 {
            return Err(FieldSeparationIntervalError::EmptyInterval {
                lower: lower.0,
                upper: upper.0,
            });
        }
        Ok(Self {
            source_id,
            reference_group_id,
            target_group_id,
            lower: AffineBoundSide {
                bound: canonical_zero(lower.0),
                configuration: lower.1,
            },
            upper: AffineBoundSide {
                bound: canonical_zero(upper.0),
                configuration: upper.1,
            },
        })
    }

    /// Returns the stable caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the ordered reference shared level set.
    pub fn reference_group_id(&self) -> &GroupId {
        &self.reference_group_id
    }

    /// Returns the ordered target shared level set.
    pub fn target_group_id(&self) -> &GroupId {
        &self.target_group_id
    }

    /// Returns the finite lower bound on `target - reference`.
    pub fn lower_bound(&self) -> f64 {
        self.lower.bound
    }

    /// Returns the finite upper bound on `target - reference`.
    pub fn upper_bound(&self) -> f64 {
        self.upper.bound
    }

    /// Reports whether either interval side owns a violation channel.
    pub fn is_soft(&self) -> bool {
        self.lower.configuration.is_soft() || self.upper.configuration.is_soft()
    }

    pub(crate) fn lower(&self) -> &AffineBoundSide {
        &self.lower
    }

    pub(crate) fn upper(&self) -> &AffineBoundSide {
        &self.upper
    }
}

/// A rejected Field Separation Interval.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FieldSeparationIntervalError {
    /// Reference and target named the same shared level set.
    SelfReference { group_id: GroupId },
    /// A lower or upper field-value bound was NaN or infinite.
    NonFiniteBound,
    /// The closed interval had its lower endpoint above its upper endpoint.
    EmptyInterval { lower: f64, upper: f64 },
}

impl fmt::Display for FieldSeparationIntervalError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SelfReference { group_id } => write!(
                formatter,
                "field-separation reference and target both name GroupId `{group_id}`"
            ),
            Self::NonFiniteBound => formatter.write_str("field-separation bound is not finite"),
            Self::EmptyInterval { lower, upper } => write!(
                formatter,
                "field-separation interval [{lower}, {upper}] is empty"
            ),
        }
    }
}

impl Error for FieldSeparationIntervalError {}

/// The explicitly declared field-value side of a shared level set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PointToLevelSetSide {
    /// The point value is greater than the shared level value.
    Increasing,
    /// The point value is less than the shared level value.
    Decreasing,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PointToLevelSetOrientation {
    lower_bounded: bool,
    bound_multiplier: f64,
    recovered_offset_multiplier: f64,
    semantic_role: &'static str,
}

impl PointToLevelSetOrientation {
    pub(crate) fn is_lower_bounded(self) -> bool {
        self.lower_bounded
    }

    pub(crate) fn bound(self, minimum_offset: MinimumFieldOffset) -> f64 {
        self.bound_multiplier * minimum_offset.value()
    }

    pub(crate) fn recovered_field_offset(self, signed_difference: f64) -> f64 {
        self.recovered_offset_multiplier * signed_difference
    }

    pub(crate) fn semantic_role(self) -> SemanticRolePath {
        SemanticRolePath::new(self.semantic_role)
    }
}

/// A finite, strictly positive point-to-level field-value difference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MinimumFieldOffset {
    value: f64,
}

impl MinimumFieldOffset {
    /// Creates a checked minimum field offset.
    pub fn try_new(value: f64) -> Result<Self, MinimumFieldOffsetError> {
        if !value.is_finite() {
            return Err(MinimumFieldOffsetError::NotFinite);
        }
        if value <= 0.0 {
            return Err(MinimumFieldOffsetError::NotPositive);
        }
        Ok(Self { value })
    }

    /// Returns the offset in the problem's field-value units.
    pub fn value(self) -> f64 {
        self.value
    }
}

/// A rejected minimum field offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MinimumFieldOffsetError {
    /// The supplied value was NaN or infinite.
    NotFinite,
    /// The supplied value was zero or negative.
    NotPositive,
}

impl fmt::Display for MinimumFieldOffsetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => formatter.write_str("minimum field offset is not finite"),
            Self::NotPositive => formatter.write_str("minimum field offset is not positive"),
        }
    }
}

impl Error for MinimumFieldOffsetError {}

/// A finite sampled point constrained to one explicit side of a shared level set.
///
/// This relation is guaranteed only at its stored support. Its offset has
/// field-value units and is not a distance or physical thickness.
#[derive(Debug, Clone, PartialEq)]
pub struct PointToLevelSetRelation {
    source_id: SourceId,
    location: Point3,
    group_id: GroupId,
    side: PointToLevelSetSide,
    minimum_offset: MinimumFieldOffset,
    configuration: AffineBoundConfiguration,
}

impl PointToLevelSetRelation {
    /// Creates one hard point-to-level-set side relation.
    pub fn hard(
        source_id: SourceId,
        location: Point3,
        group_id: GroupId,
        side: PointToLevelSetSide,
        minimum_offset: MinimumFieldOffset,
    ) -> Self {
        Self::new(
            source_id,
            location,
            group_id,
            side,
            minimum_offset,
            AffineBoundConfiguration::Hard,
        )
    }

    /// Creates a soft point-side relation with a quadratic violation loss.
    pub fn with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        group_id: GroupId,
        side: PointToLevelSetSide,
        minimum_offset: MinimumFieldOffset,
        penalty: QuadraticPenalty,
    ) -> Self {
        Self::new(
            source_id,
            location,
            group_id,
            side,
            minimum_offset,
            AffineBoundConfiguration::QuadraticPenalty(penalty),
        )
    }

    /// Creates a soft point-side relation with a linear violation loss.
    pub fn with_linear_violation_penalty(
        source_id: SourceId,
        location: Point3,
        group_id: GroupId,
        side: PointToLevelSetSide,
        minimum_offset: MinimumFieldOffset,
        penalty: LinearViolationPenalty,
    ) -> Self {
        Self::new(
            source_id,
            location,
            group_id,
            side,
            minimum_offset,
            AffineBoundConfiguration::LinearViolationPenalty(penalty),
        )
    }

    fn new(
        source_id: SourceId,
        location: Point3,
        group_id: GroupId,
        side: PointToLevelSetSide,
        minimum_offset: MinimumFieldOffset,
        configuration: AffineBoundConfiguration,
    ) -> Self {
        Self {
            source_id,
            location,
            group_id,
            side,
            minimum_offset,
            configuration,
        }
    }

    /// Returns the stable caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the finite sampled location in the declared input frame.
    pub fn location(&self) -> Point3 {
        self.location
    }

    /// Returns the referenced shared level set.
    pub fn group_id(&self) -> &GroupId {
        &self.group_id
    }

    /// Returns the explicitly declared field-value side.
    pub fn side(&self) -> PointToLevelSetSide {
        self.side
    }

    pub(crate) fn orientation(&self) -> PointToLevelSetOrientation {
        match self.side {
            PointToLevelSetSide::Increasing => PointToLevelSetOrientation {
                lower_bounded: true,
                bound_multiplier: 1.0,
                recovered_offset_multiplier: 1.0,
                semantic_role: "point-to-level-set/increasing/minimum-field-offset",
            },
            PointToLevelSetSide::Decreasing => PointToLevelSetOrientation {
                lower_bounded: false,
                bound_multiplier: -1.0,
                recovered_offset_multiplier: -1.0,
                semantic_role: "point-to-level-set/decreasing/minimum-field-offset",
            },
        }
    }

    /// Returns the strictly positive field-value offset.
    pub fn minimum_offset(&self) -> MinimumFieldOffset {
        self.minimum_offset
    }

    /// Reports whether this relation owns an explicit violation channel.
    pub fn is_soft(&self) -> bool {
        self.configuration.is_soft()
    }

    pub(crate) fn configuration(&self) -> AffineBoundConfiguration {
        self.configuration
    }
}

/// The explicit mapping from stratigraphic age to scalar-field value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum StratigraphicFieldDirection {
    /// Field value increases toward geologically younger horizons.
    TowardYounger,
    /// Field value increases toward geologically older horizons.
    TowardOlder,
}

/// A finite, strictly positive shared-level-set difference in field-value units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MinimumFieldSeparation {
    value: f64,
}

impl MinimumFieldSeparation {
    /// Creates a checked minimum field separation.
    pub fn try_new(value: f64) -> Result<Self, MinimumFieldSeparationError> {
        if !value.is_finite() {
            return Err(MinimumFieldSeparationError::NotFinite);
        }
        if value <= 0.0 {
            return Err(MinimumFieldSeparationError::NotPositive);
        }
        Ok(Self { value })
    }

    /// Returns the separation in the problem's field-value units.
    pub fn value(self) -> f64 {
        self.value
    }
}

/// A rejected minimum field separation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MinimumFieldSeparationError {
    /// The supplied value was NaN or infinite.
    NotFinite,
    /// The supplied value was zero or negative.
    NotPositive,
}

impl fmt::Display for MinimumFieldSeparationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => formatter.write_str("minimum field separation is not finite"),
            Self::NotPositive => formatter.write_str("minimum field separation is not positive"),
        }
    }
}

impl Error for MinimumFieldSeparationError {}

/// The caller-declared semantic kind of a relation between shared level sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SharedLevelSetRelationKind {
    /// One named level is geologically younger than another.
    YoungerThan,
    /// One named level is geologically older than another.
    OlderThan,
    /// One field level is non-strictly no greater than another.
    FieldLevelOrder,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum SharedLevelSetRelationInput {
    YoungerThan {
        source_id: SourceId,
        younger_group_id: GroupId,
        older_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        configuration: AffineBoundConfiguration,
    },
    OlderThan {
        source_id: SourceId,
        older_group_id: GroupId,
        younger_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        configuration: AffineBoundConfiguration,
    },
    FieldLevelOrder {
        source_id: SourceId,
        lower_group_id: GroupId,
        upper_group_id: GroupId,
        configuration: AffineBoundConfiguration,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SharedLevelSetRelationOrientation {
    pub(crate) lower_field_group_id: GroupId,
    pub(crate) upper_field_group_id: GroupId,
    pub(crate) field_direction: Option<StratigraphicFieldDirection>,
    pub(crate) minimum_separation: Option<MinimumFieldSeparation>,
}

impl SharedLevelSetRelationOrientation {
    pub(crate) fn required_difference(&self) -> f64 {
        self.minimum_separation
            .map(MinimumFieldSeparation::value)
            .unwrap_or(0.0)
    }
}

impl SharedLevelSetRelationInput {
    fn younger_than(
        source_id: SourceId,
        younger_group_id: GroupId,
        older_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        configuration: AffineBoundConfiguration,
    ) -> Self {
        Self::YoungerThan {
            source_id,
            younger_group_id,
            older_group_id,
            minimum_separation,
            configuration,
        }
    }

    fn older_than(
        source_id: SourceId,
        older_group_id: GroupId,
        younger_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        configuration: AffineBoundConfiguration,
    ) -> Self {
        Self::OlderThan {
            source_id,
            older_group_id,
            younger_group_id,
            minimum_separation,
            configuration,
        }
    }

    fn field_level_order(
        source_id: SourceId,
        lower_group_id: GroupId,
        upper_group_id: GroupId,
        configuration: AffineBoundConfiguration,
    ) -> Self {
        Self::FieldLevelOrder {
            source_id,
            lower_group_id,
            upper_group_id,
            configuration,
        }
    }

    pub(crate) fn source_id(&self) -> &SourceId {
        match self {
            Self::YoungerThan { source_id, .. }
            | Self::OlderThan { source_id, .. }
            | Self::FieldLevelOrder { source_id, .. } => source_id,
        }
    }

    pub(crate) fn declared_group_ids(&self) -> [&GroupId; 2] {
        match self {
            Self::YoungerThan {
                younger_group_id,
                older_group_id,
                ..
            } => [younger_group_id, older_group_id],
            Self::OlderThan {
                older_group_id,
                younger_group_id,
                ..
            } => [older_group_id, younger_group_id],
            Self::FieldLevelOrder {
                lower_group_id,
                upper_group_id,
                ..
            } => [lower_group_id, upper_group_id],
        }
    }

    pub(crate) fn kind(&self) -> SharedLevelSetRelationKind {
        match self {
            Self::YoungerThan { .. } => SharedLevelSetRelationKind::YoungerThan,
            Self::OlderThan { .. } => SharedLevelSetRelationKind::OlderThan,
            Self::FieldLevelOrder { .. } => SharedLevelSetRelationKind::FieldLevelOrder,
        }
    }

    pub(crate) fn younger_group_id(&self) -> Option<&GroupId> {
        match self {
            Self::YoungerThan {
                younger_group_id, ..
            }
            | Self::OlderThan {
                younger_group_id, ..
            } => Some(younger_group_id),
            Self::FieldLevelOrder { .. } => None,
        }
    }

    pub(crate) fn older_group_id(&self) -> Option<&GroupId> {
        match self {
            Self::YoungerThan { older_group_id, .. } | Self::OlderThan { older_group_id, .. } => {
                Some(older_group_id)
            }
            Self::FieldLevelOrder { .. } => None,
        }
    }

    pub(crate) fn lower_group_id(&self) -> Option<&GroupId> {
        match self {
            Self::FieldLevelOrder { lower_group_id, .. } => Some(lower_group_id),
            Self::YoungerThan { .. } | Self::OlderThan { .. } => None,
        }
    }

    pub(crate) fn upper_group_id(&self) -> Option<&GroupId> {
        match self {
            Self::FieldLevelOrder { upper_group_id, .. } => Some(upper_group_id),
            Self::YoungerThan { .. } | Self::OlderThan { .. } => None,
        }
    }

    pub(crate) fn minimum_separation(&self) -> Option<MinimumFieldSeparation> {
        match self {
            Self::YoungerThan {
                minimum_separation, ..
            }
            | Self::OlderThan {
                minimum_separation, ..
            } => Some(*minimum_separation),
            Self::FieldLevelOrder { .. } => None,
        }
    }

    pub(crate) fn configuration(&self) -> AffineBoundConfiguration {
        match self {
            Self::YoungerThan { configuration, .. }
            | Self::OlderThan { configuration, .. }
            | Self::FieldLevelOrder { configuration, .. } => *configuration,
        }
    }

    pub(crate) fn is_soft(&self) -> bool {
        self.configuration().is_soft()
    }

    pub(crate) fn is_stratigraphic_age_relation(&self) -> bool {
        matches!(self, Self::YoungerThan { .. } | Self::OlderThan { .. })
    }

    pub(crate) fn is_field_level_order(&self) -> bool {
        matches!(self, Self::FieldLevelOrder { .. })
    }

    pub(crate) fn semantic_role(&self) -> SemanticRolePath {
        match self {
            Self::YoungerThan { .. } => {
                SemanticRolePath::new("younger-than/minimum-field-separation")
            }
            Self::OlderThan { .. } => SemanticRolePath::new("older-than/minimum-field-separation"),
            Self::FieldLevelOrder { .. } => SemanticRolePath::new("field-level-order/non-strict"),
        }
    }

    pub(crate) fn orientation(
        &self,
        direction: Option<StratigraphicFieldDirection>,
    ) -> SharedLevelSetRelationOrientation {
        let age_direction =
            || direction.expect("build requires a field direction for every age relation");
        match self {
            Self::FieldLevelOrder {
                lower_group_id,
                upper_group_id,
                ..
            } => SharedLevelSetRelationOrientation {
                lower_field_group_id: lower_group_id.clone(),
                upper_field_group_id: upper_group_id.clone(),
                field_direction: None,
                minimum_separation: None,
            },
            Self::YoungerThan {
                younger_group_id,
                older_group_id,
                minimum_separation,
                ..
            } => {
                let direction = age_direction();
                let (lower, upper) = match direction {
                    StratigraphicFieldDirection::TowardYounger => {
                        (older_group_id, younger_group_id)
                    }
                    StratigraphicFieldDirection::TowardOlder => (younger_group_id, older_group_id),
                };
                SharedLevelSetRelationOrientation {
                    lower_field_group_id: lower.clone(),
                    upper_field_group_id: upper.clone(),
                    field_direction: Some(direction),
                    minimum_separation: Some(*minimum_separation),
                }
            }
            Self::OlderThan {
                older_group_id,
                younger_group_id,
                minimum_separation,
                ..
            } => {
                let direction = age_direction();
                let (lower, upper) = match direction {
                    StratigraphicFieldDirection::TowardYounger => {
                        (older_group_id, younger_group_id)
                    }
                    StratigraphicFieldDirection::TowardOlder => (younger_group_id, older_group_id),
                };
                SharedLevelSetRelationOrientation {
                    lower_field_group_id: lower.clone(),
                    upper_field_group_id: upper.clone(),
                    field_direction: Some(direction),
                    minimum_separation: Some(*minimum_separation),
                }
            }
        }
    }
}

/// A strict stratigraphic statement that one horizon is younger than another.
#[derive(Debug, Clone, PartialEq)]
pub struct YoungerThan(SharedLevelSetRelationInput);

impl YoungerThan {
    /// Creates one hard younger-than relation.
    pub fn hard(
        source_id: SourceId,
        younger_group_id: GroupId,
        older_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
    ) -> Self {
        Self(SharedLevelSetRelationInput::younger_than(
            source_id,
            younger_group_id,
            older_group_id,
            minimum_separation,
            AffineBoundConfiguration::Hard,
        ))
    }

    /// Creates a soft younger-than relation with a quadratic violation loss.
    pub fn with_quadratic_penalty(
        source_id: SourceId,
        younger_group_id: GroupId,
        older_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        penalty: QuadraticPenalty,
    ) -> Self {
        Self(SharedLevelSetRelationInput::younger_than(
            source_id,
            younger_group_id,
            older_group_id,
            minimum_separation,
            AffineBoundConfiguration::QuadraticPenalty(penalty),
        ))
    }

    /// Creates a soft younger-than relation with a linear violation loss.
    pub fn with_linear_violation_penalty(
        source_id: SourceId,
        younger_group_id: GroupId,
        older_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        penalty: LinearViolationPenalty,
    ) -> Self {
        Self(SharedLevelSetRelationInput::younger_than(
            source_id,
            younger_group_id,
            older_group_id,
            minimum_separation,
            AffineBoundConfiguration::LinearViolationPenalty(penalty),
        ))
    }

    /// Returns the stable caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        self.0.source_id()
    }

    /// Returns the geologically younger shared level set.
    pub fn younger_group_id(&self) -> &GroupId {
        self.0
            .younger_group_id()
            .expect("YoungerThan always owns a younger group")
    }

    /// Returns the geologically older shared level set.
    pub fn older_group_id(&self) -> &GroupId {
        self.0
            .older_group_id()
            .expect("YoungerThan always owns an older group")
    }

    /// Returns the required strict field-value difference.
    pub fn minimum_separation(&self) -> MinimumFieldSeparation {
        self.0
            .minimum_separation()
            .expect("an age relation always owns a separation")
    }

    /// Reports whether this relation owns an explicit violation channel.
    pub fn is_soft(&self) -> bool {
        self.0.is_soft()
    }
}

/// A strict stratigraphic statement that one horizon is older than another.
#[derive(Debug, Clone, PartialEq)]
pub struct OlderThan(SharedLevelSetRelationInput);

impl OlderThan {
    /// Creates one hard older-than relation.
    pub fn hard(
        source_id: SourceId,
        older_group_id: GroupId,
        younger_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
    ) -> Self {
        Self(SharedLevelSetRelationInput::older_than(
            source_id,
            older_group_id,
            younger_group_id,
            minimum_separation,
            AffineBoundConfiguration::Hard,
        ))
    }

    /// Creates a soft older-than relation with a quadratic violation loss.
    pub fn with_quadratic_penalty(
        source_id: SourceId,
        older_group_id: GroupId,
        younger_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        penalty: QuadraticPenalty,
    ) -> Self {
        Self(SharedLevelSetRelationInput::older_than(
            source_id,
            older_group_id,
            younger_group_id,
            minimum_separation,
            AffineBoundConfiguration::QuadraticPenalty(penalty),
        ))
    }

    /// Creates a soft older-than relation with a linear violation loss.
    pub fn with_linear_violation_penalty(
        source_id: SourceId,
        older_group_id: GroupId,
        younger_group_id: GroupId,
        minimum_separation: MinimumFieldSeparation,
        penalty: LinearViolationPenalty,
    ) -> Self {
        Self(SharedLevelSetRelationInput::older_than(
            source_id,
            older_group_id,
            younger_group_id,
            minimum_separation,
            AffineBoundConfiguration::LinearViolationPenalty(penalty),
        ))
    }

    /// Returns the stable caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        self.0.source_id()
    }

    /// Returns the geologically older shared level set.
    pub fn older_group_id(&self) -> &GroupId {
        self.0
            .older_group_id()
            .expect("OlderThan always owns an older group")
    }

    /// Returns the geologically younger shared level set.
    pub fn younger_group_id(&self) -> &GroupId {
        self.0
            .younger_group_id()
            .expect("OlderThan always owns a younger group")
    }

    /// Returns the required strict field-value difference.
    pub fn minimum_separation(&self) -> MinimumFieldSeparation {
        self.0
            .minimum_separation()
            .expect("an age relation always owns a separation")
    }

    /// Reports whether this relation owns an explicit violation channel.
    pub fn is_soft(&self) -> bool {
        self.0.is_soft()
    }
}

/// A non-strict direct ordering of two shared scalar-field levels.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldLevelOrder(SharedLevelSetRelationInput);

impl FieldLevelOrder {
    /// Creates the hard relation `lower_group <= upper_group`.
    pub fn hard(source_id: SourceId, lower_group_id: GroupId, upper_group_id: GroupId) -> Self {
        Self(SharedLevelSetRelationInput::field_level_order(
            source_id,
            lower_group_id,
            upper_group_id,
            AffineBoundConfiguration::Hard,
        ))
    }

    /// Creates a soft non-strict order with a quadratic violation loss.
    pub fn with_quadratic_penalty(
        source_id: SourceId,
        lower_group_id: GroupId,
        upper_group_id: GroupId,
        penalty: QuadraticPenalty,
    ) -> Self {
        Self(SharedLevelSetRelationInput::field_level_order(
            source_id,
            lower_group_id,
            upper_group_id,
            AffineBoundConfiguration::QuadraticPenalty(penalty),
        ))
    }

    /// Creates a soft non-strict order with a linear violation loss.
    pub fn with_linear_violation_penalty(
        source_id: SourceId,
        lower_group_id: GroupId,
        upper_group_id: GroupId,
        penalty: LinearViolationPenalty,
    ) -> Self {
        Self(SharedLevelSetRelationInput::field_level_order(
            source_id,
            lower_group_id,
            upper_group_id,
            AffineBoundConfiguration::LinearViolationPenalty(penalty),
        ))
    }

    /// Returns the stable caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        self.0.source_id()
    }

    /// Returns the shared level set constrained to have the lower-or-equal value.
    pub fn lower_group_id(&self) -> &GroupId {
        self.0
            .lower_group_id()
            .expect("FieldLevelOrder always owns a lower group")
    }

    /// Returns the shared level set constrained to have the upper-or-equal value.
    pub fn upper_group_id(&self) -> &GroupId {
        self.0
            .upper_group_id()
            .expect("FieldLevelOrder always owns an upper group")
    }

    /// Reports whether this relation owns an explicit violation channel.
    pub fn is_soft(&self) -> bool {
        self.0.is_soft()
    }
}

impl From<QuadraticPenalty> for DirectionalDerivativeViolationPenalty {
    fn from(penalty: QuadraticPenalty) -> Self {
        Self::Quadratic(penalty)
    }
}

impl From<LinearViolationPenalty> for DirectionalDerivativeViolationPenalty {
    fn from(penalty: LinearViolationPenalty) -> Self {
        Self::Linear(penalty)
    }
}

impl From<QuadraticPenalty> for FieldValueViolationPenalty {
    fn from(penalty: QuadraticPenalty) -> Self {
        Self::Quadratic(penalty)
    }
}

impl From<LinearViolationPenalty> for FieldValueViolationPenalty {
    fn from(penalty: LinearViolationPenalty) -> Self {
        Self::Linear(penalty)
    }
}

/// A checked lower, upper, or interval bound on the field at one finite location.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldValueBound {
    source_id: SourceId,
    location: Point3,
    lower: Option<AffineBoundSide>,
    upper: Option<AffineBoundSide>,
}

impl FieldValueBound {
    /// Creates one hard finite lower bound.
    pub fn try_lower(
        source_id: SourceId,
        location: Point3,
        lower: f64,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            Some((lower, AffineBoundConfiguration::Hard)),
            None,
        )
    }

    /// Creates one hard finite upper bound.
    pub fn try_upper(
        source_id: SourceId,
        location: Point3,
        upper: f64,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            None,
            Some((upper, AffineBoundConfiguration::Hard)),
        )
    }

    /// Creates one hard closed finite interval.
    pub fn try_interval(
        source_id: SourceId,
        location: Point3,
        lower: f64,
        upper: f64,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            Some((lower, AffineBoundConfiguration::Hard)),
            Some((upper, AffineBoundConfiguration::Hard)),
        )
    }

    /// Creates one soft lower bound with a quadratic violation penalty.
    pub fn try_lower_with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        lower: f64,
        penalty: QuadraticPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            Some((lower, AffineBoundConfiguration::QuadraticPenalty(penalty))),
            None,
        )
    }

    /// Creates one soft upper bound with a quadratic violation penalty.
    pub fn try_upper_with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        upper: f64,
        penalty: QuadraticPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            None,
            Some((upper, AffineBoundConfiguration::QuadraticPenalty(penalty))),
        )
    }

    /// Creates a soft interval with independent lower and upper quadratic penalties.
    pub fn try_interval_with_quadratic_penalties(
        source_id: SourceId,
        location: Point3,
        lower: f64,
        lower_penalty: QuadraticPenalty,
        upper: f64,
        upper_penalty: QuadraticPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            Some((
                lower,
                AffineBoundConfiguration::QuadraticPenalty(lower_penalty),
            )),
            Some((
                upper,
                AffineBoundConfiguration::QuadraticPenalty(upper_penalty),
            )),
        )
    }

    /// Creates one soft lower bound with a linear violation penalty.
    pub fn try_lower_with_linear_violation_penalty(
        source_id: SourceId,
        location: Point3,
        lower: f64,
        penalty: LinearViolationPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            Some((
                lower,
                AffineBoundConfiguration::LinearViolationPenalty(penalty),
            )),
            None,
        )
    }

    /// Creates one soft upper bound with a linear violation penalty.
    pub fn try_upper_with_linear_violation_penalty(
        source_id: SourceId,
        location: Point3,
        upper: f64,
        penalty: LinearViolationPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            None,
            Some((
                upper,
                AffineBoundConfiguration::LinearViolationPenalty(penalty),
            )),
        )
    }

    /// Creates a soft interval with independent lower and upper linear penalties.
    pub fn try_interval_with_linear_violation_penalties(
        source_id: SourceId,
        location: Point3,
        lower: f64,
        lower_penalty: LinearViolationPenalty,
        upper: f64,
        upper_penalty: LinearViolationPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        Self::new(
            source_id,
            location,
            Some((
                lower,
                AffineBoundConfiguration::LinearViolationPenalty(lower_penalty),
            )),
            Some((
                upper,
                AffineBoundConfiguration::LinearViolationPenalty(upper_penalty),
            )),
        )
    }

    /// Creates a soft interval whose two sides independently select a legal loss.
    pub fn try_interval_with_violation_penalties(
        source_id: SourceId,
        location: Point3,
        lower: f64,
        lower_penalty: FieldValueViolationPenalty,
        upper: f64,
        upper_penalty: FieldValueViolationPenalty,
    ) -> Result<Self, FieldValueBoundError> {
        let configuration = |penalty| match penalty {
            FieldValueViolationPenalty::Quadratic(penalty) => {
                AffineBoundConfiguration::QuadraticPenalty(penalty)
            }
            FieldValueViolationPenalty::Linear(penalty) => {
                AffineBoundConfiguration::LinearViolationPenalty(penalty)
            }
        };
        Self::new(
            source_id,
            location,
            Some((lower, configuration(lower_penalty))),
            Some((upper, configuration(upper_penalty))),
        )
    }

    fn new(
        source_id: SourceId,
        location: Point3,
        lower: Option<(f64, AffineBoundConfiguration)>,
        upper: Option<(f64, AffineBoundConfiguration)>,
    ) -> Result<Self, FieldValueBoundError> {
        if lower
            .iter()
            .chain(&upper)
            .any(|(bound, _)| !bound.is_finite())
        {
            return Err(FieldValueBoundError::NonFiniteBound);
        }
        if let (Some((lower, _)), Some((upper, _))) = (lower, upper) {
            if lower > upper {
                return Err(FieldValueBoundError::EmptyInterval { lower, upper });
            }
        }
        let side = |(bound, configuration): (f64, AffineBoundConfiguration)| AffineBoundSide {
            bound: canonical_zero(bound),
            configuration,
        };
        Ok(Self {
            source_id,
            location,
            lower: lower.map(side),
            upper: upper.map(side),
        })
    }

    /// Returns the stable caller-owned identity of this relation.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the finite support location in the declared input frame.
    pub fn location(&self) -> Point3 {
        self.location
    }

    /// Returns the lower field-value bound when present.
    pub fn lower_bound(&self) -> Option<f64> {
        self.lower.as_ref().map(|side| side.bound)
    }

    /// Returns the upper field-value bound when present.
    pub fn upper_bound(&self) -> Option<f64> {
        self.upper.as_ref().map(|side| side.bound)
    }

    /// Reports whether this relation owns explicit violation channels.
    pub fn is_soft(&self) -> bool {
        self.lower
            .iter()
            .chain(&self.upper)
            .any(|side| side.configuration.is_soft())
    }

    pub(crate) fn lower(&self) -> Option<&AffineBoundSide> {
        self.lower.as_ref()
    }

    pub(crate) fn upper(&self) -> Option<&AffineBoundSide> {
        self.upper.as_ref()
    }
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

/// A rejected Field Value Bound.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FieldValueBoundError {
    /// A lower or upper bound was NaN or infinite.
    NonFiniteBound,
    /// A closed interval had its lower endpoint above its upper endpoint.
    EmptyInterval { lower: f64, upper: f64 },
}

impl fmt::Display for FieldValueBoundError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteBound => formatter.write_str("field-value bound is not finite"),
            Self::EmptyInterval { lower, upper } => {
                write!(
                    formatter,
                    "field-value interval [{lower}, {upper}] is empty"
                )
            }
        }
    }
}

impl Error for FieldValueBoundError {}

/// A checked lower, upper, or interval bound on one directional derivative.
///
/// The derivative is taken along the stored oriented unit direction in the
/// problem's physical input coordinates. Bounds therefore use field-value per
/// length units; they are not angular or numerical tolerances and do not imply
/// a complete gradient magnitude.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectionalDerivativeInterval {
    source_id: SourceId,
    location: Point3,
    direction: Vector3,
    lower: Option<AffineBoundSide>,
    upper: Option<AffineBoundSide>,
}

impl DirectionalDerivativeInterval {
    /// Creates one hard finite lower directional-derivative bound.
    pub fn try_lower(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            Some((lower, AffineBoundConfiguration::Hard)),
            None,
        )
    }

    /// Creates one hard finite upper directional-derivative bound.
    pub fn try_upper(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        upper: f64,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            None,
            Some((upper, AffineBoundConfiguration::Hard)),
        )
    }

    /// Creates one hard closed finite directional-derivative interval.
    pub fn try_interval(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
        upper: f64,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            Some((lower, AffineBoundConfiguration::Hard)),
            Some((upper, AffineBoundConfiguration::Hard)),
        )
    }

    /// Creates one soft lower bound with a quadratic violation penalty.
    pub fn try_lower_with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
        penalty: QuadraticPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            Some((lower, AffineBoundConfiguration::QuadraticPenalty(penalty))),
            None,
        )
    }

    /// Creates one soft upper bound with a quadratic violation penalty.
    pub fn try_upper_with_quadratic_penalty(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        upper: f64,
        penalty: QuadraticPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            None,
            Some((upper, AffineBoundConfiguration::QuadraticPenalty(penalty))),
        )
    }

    /// Creates a soft interval with independent quadratic penalties per side.
    #[allow(clippy::too_many_arguments)]
    pub fn try_interval_with_quadratic_penalties(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
        lower_penalty: QuadraticPenalty,
        upper: f64,
        upper_penalty: QuadraticPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            Some((
                lower,
                AffineBoundConfiguration::QuadraticPenalty(lower_penalty),
            )),
            Some((
                upper,
                AffineBoundConfiguration::QuadraticPenalty(upper_penalty),
            )),
        )
    }

    /// Creates one soft lower bound with a linear violation penalty.
    pub fn try_lower_with_linear_violation_penalty(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
        penalty: LinearViolationPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            Some((
                lower,
                AffineBoundConfiguration::LinearViolationPenalty(penalty),
            )),
            None,
        )
    }

    /// Creates one soft upper bound with a linear violation penalty.
    pub fn try_upper_with_linear_violation_penalty(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        upper: f64,
        penalty: LinearViolationPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            None,
            Some((
                upper,
                AffineBoundConfiguration::LinearViolationPenalty(penalty),
            )),
        )
    }

    /// Creates a soft interval with independent linear penalties per side.
    #[allow(clippy::too_many_arguments)]
    pub fn try_interval_with_linear_violation_penalties(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
        lower_penalty: LinearViolationPenalty,
        upper: f64,
        upper_penalty: LinearViolationPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        Self::new(
            source_id,
            location,
            direction,
            Some((
                lower,
                AffineBoundConfiguration::LinearViolationPenalty(lower_penalty),
            )),
            Some((
                upper,
                AffineBoundConfiguration::LinearViolationPenalty(upper_penalty),
            )),
        )
    }

    /// Creates a soft interval whose sides independently select a legal loss.
    #[allow(clippy::too_many_arguments)]
    pub fn try_interval_with_violation_penalties(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: f64,
        lower_penalty: DirectionalDerivativeViolationPenalty,
        upper: f64,
        upper_penalty: DirectionalDerivativeViolationPenalty,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        let configuration = |penalty| match penalty {
            DirectionalDerivativeViolationPenalty::Quadratic(penalty) => {
                AffineBoundConfiguration::QuadraticPenalty(penalty)
            }
            DirectionalDerivativeViolationPenalty::Linear(penalty) => {
                AffineBoundConfiguration::LinearViolationPenalty(penalty)
            }
        };
        Self::new(
            source_id,
            location,
            direction,
            Some((lower, configuration(lower_penalty))),
            Some((upper, configuration(upper_penalty))),
        )
    }

    fn new(
        source_id: SourceId,
        location: Point3,
        direction: Vector3,
        lower: Option<(f64, AffineBoundConfiguration)>,
        upper: Option<(f64, AffineBoundConfiguration)>,
    ) -> Result<Self, DirectionalDerivativeIntervalError> {
        if lower
            .iter()
            .chain(&upper)
            .any(|(bound, _)| !bound.is_finite())
        {
            return Err(DirectionalDerivativeIntervalError::NonFiniteBound);
        }
        if let (Some((lower, _)), Some((upper, _))) = (lower, upper) {
            if lower > upper {
                return Err(DirectionalDerivativeIntervalError::EmptyInterval { lower, upper });
            }
        }
        let direction = normalize_direction(direction)
            .ok_or(DirectionalDerivativeIntervalError::ZeroDirection)?;
        let side = |(bound, configuration)| AffineBoundSide {
            bound: canonical_zero(bound),
            configuration,
        };
        Ok(Self {
            source_id,
            location,
            direction,
            lower: lower.map(side),
            upper: upper.map(side),
        })
    }

    /// Returns the stable caller-owned relation identity.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the finite support location in the declared input frame.
    pub fn location(&self) -> Point3 {
        self.location
    }

    /// Returns the oriented unit direction in physical input coordinates.
    pub fn direction(&self) -> Vector3 {
        self.direction
    }

    /// Returns the lower directional-derivative bound when present.
    pub fn lower_bound(&self) -> Option<f64> {
        self.lower.as_ref().map(|side| side.bound)
    }

    /// Returns the upper directional-derivative bound when present.
    pub fn upper_bound(&self) -> Option<f64> {
        self.upper.as_ref().map(|side| side.bound)
    }

    /// Reports whether this relation owns explicit violation channels.
    pub fn is_soft(&self) -> bool {
        self.lower
            .iter()
            .chain(&self.upper)
            .any(|side| side.configuration.is_soft())
    }

    pub(crate) fn lower(&self) -> Option<&AffineBoundSide> {
        self.lower.as_ref()
    }

    pub(crate) fn upper(&self) -> Option<&AffineBoundSide> {
        self.upper.as_ref()
    }
}

/// A rejected Directional Derivative Interval.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum DirectionalDerivativeIntervalError {
    /// The physical direction was the zero vector.
    ZeroDirection,
    /// A lower or upper derivative bound was NaN or infinite.
    NonFiniteBound,
    /// A closed interval had its lower endpoint above its upper endpoint.
    EmptyInterval { lower: f64, upper: f64 },
}

impl fmt::Display for DirectionalDerivativeIntervalError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroDirection => formatter.write_str("directional-derivative direction is zero"),
            Self::NonFiniteBound => {
                formatter.write_str("directional-derivative bound is not finite")
            }
            Self::EmptyInterval { lower, upper } => write!(
                formatter,
                "directional-derivative interval [{lower}, {upper}] is empty"
            ),
        }
    }
}

impl Error for DirectionalDerivativeIntervalError {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AffineBoundSide {
    pub(crate) bound: f64,
    pub(crate) configuration: AffineBoundConfiguration,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum AffineBoundConfiguration {
    Hard,
    QuadraticPenalty(QuadraticPenalty),
    LinearViolationPenalty(LinearViolationPenalty),
}

impl AffineBoundConfiguration {
    pub(crate) fn is_soft(self) -> bool {
        !matches!(self, Self::Hard)
    }
}

/// One immutable member of a shared field-value group.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedLevelSetMember {
    source_id: SourceId,
    location: Point3,
}

impl SharedLevelSetMember {
    /// Returns the caller-owned identity of this member relation.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the member location in the declared input frame.
    pub fn location(&self) -> Point3 {
        self.location
    }
}

#[derive(Debug)]
struct GroupDraft {
    group_id: GroupId,
    members: Vec<SharedLevelSetMember>,
    source_ids: BTreeSet<SourceId>,
}

impl GroupDraft {
    fn new(group_id: GroupId) -> Self {
        Self {
            group_id,
            members: Vec::new(),
            source_ids: BTreeSet::new(),
        }
    }

    fn add_member(
        &mut self,
        source_id: SourceId,
        location: Point3,
    ) -> Result<(), GroupMemberAddError> {
        if self.source_ids.contains(&source_id) {
            return Err(GroupMemberAddError::DuplicateSourceId { source_id });
        }
        self.source_ids.insert(source_id.clone());
        self.members.push(SharedLevelSetMember {
            source_id,
            location,
        });
        Ok(())
    }

    fn finish(mut self) -> Result<CompletedGroup, GroupBuildError> {
        if self.members.is_empty() {
            return Err(GroupBuildError::EmptyGroup);
        }
        self.members
            .sort_by(|left, right| left.source_id.cmp(&right.source_id));
        Ok(CompletedGroup {
            group_id: self.group_id,
            members: self.members,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
struct CompletedGroup {
    group_id: GroupId,
    members: Vec<SharedLevelSetMember>,
}

/// Atomically constructs a general mathematical shared level set.
#[derive(Debug)]
pub struct SharedLevelSetBuilder {
    draft: GroupDraft,
}

impl SharedLevelSetBuilder {
    /// Starts an incomplete group with its stable identity.
    pub fn new(group_id: GroupId) -> Self {
        Self {
            draft: GroupDraft::new(group_id),
        }
    }

    /// Adds one complete hard member to this group draft.
    pub fn add_member(
        &mut self,
        source_id: SourceId,
        location: Point3,
    ) -> Result<(), GroupMemberAddError> {
        self.draft.add_member(source_id, location)
    }

    /// Finishes a non-empty immutable group.
    pub fn build(self) -> Result<SharedLevelSet, GroupBuildError> {
        self.draft.finish().map(SharedLevelSet)
    }
}

/// A general group of positions sharing one unknown semantic field value.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedLevelSet(CompletedGroup);

impl SharedLevelSet {
    /// Returns this group's caller-owned stable identity.
    pub fn group_id(&self) -> &GroupId {
        &self.0.group_id
    }

    /// Returns all members in stable SourceId order.
    pub fn members(&self) -> &[SharedLevelSetMember] {
        &self.0.members
    }
}

/// Atomically constructs a geological horizon.
#[derive(Debug)]
pub struct HorizonBuilder {
    draft: GroupDraft,
}

impl HorizonBuilder {
    /// Starts an incomplete horizon with its stable identity.
    pub fn new(group_id: GroupId) -> Self {
        Self {
            draft: GroupDraft::new(group_id),
        }
    }

    /// Adds one complete hard horizon member to this group draft.
    pub fn add_member(
        &mut self,
        source_id: SourceId,
        location: Point3,
    ) -> Result<(), GroupMemberAddError> {
        self.draft.add_member(source_id, location)
    }

    /// Finishes a non-empty immutable horizon.
    pub fn build(self) -> Result<Horizon, GroupBuildError> {
        self.draft.finish().map(Horizon)
    }
}

/// A stratigraphic shared level set with one unknown horizon field value.
#[derive(Debug, Clone, PartialEq)]
pub struct Horizon(CompletedGroup);

impl Horizon {
    /// Returns this horizon's caller-owned stable identity.
    pub fn group_id(&self) -> &GroupId {
        &self.0.group_id
    }

    /// Returns all horizon members in stable SourceId order.
    pub fn members(&self) -> &[SharedLevelSetMember] {
        &self.0.members
    }
}

/// A rejected mutation of an incomplete group draft.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum GroupMemberAddError {
    /// A member reused another member's SourceId.
    DuplicateSourceId { source_id: SourceId },
}

impl fmt::Display for GroupMemberAddError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateSourceId { source_id } => {
                write!(formatter, "duplicate group member SourceId `{source_id}`")
            }
        }
    }
}

impl Error for GroupMemberAddError {}

/// A rejected attempt to finish an incomplete group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum GroupBuildError {
    /// Shared field-value semantics require at least one member.
    EmptyGroup,
}

impl fmt::Display for GroupBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGroup => formatter.write_str("a shared level set cannot be empty"),
        }
    }
}

impl Error for GroupBuildError {}

/// Explicit convention selecting the additive constant of a scalar field.
#[derive(Debug, Clone, PartialEq)]
pub struct AdditiveFieldGauge {
    source_id: SourceId,
    reference: AdditiveFieldGaugeReference,
    value: f64,
}

impl AdditiveFieldGauge {
    /// Selects the field value at one point without claiming an observation.
    pub fn at_point(source_id: SourceId, point: Point3, value: f64) -> Result<Self, GaugeError> {
        Self::new(source_id, AdditiveFieldGaugeReference::Point(point), value)
    }

    /// Selects one declared shared level set's semantic field value.
    pub fn at_level_set(
        source_id: SourceId,
        group_id: GroupId,
        value: f64,
    ) -> Result<Self, GaugeError> {
        Self::new(
            source_id,
            AdditiveFieldGaugeReference::LevelSet(group_id),
            value,
        )
    }

    fn new(
        source_id: SourceId,
        reference: AdditiveFieldGaugeReference,
        value: f64,
    ) -> Result<Self, GaugeError> {
        if !value.is_finite() {
            return Err(GaugeError::NonFiniteFieldValue);
        }
        Ok(Self {
            source_id,
            reference,
            value,
        })
    }

    /// Returns the stable source identity of this convention.
    pub fn source_id(&self) -> &SourceId {
        &self.source_id
    }

    /// Returns the selected absolute field value.
    pub fn value(&self) -> f64 {
        self.value
    }

    pub(crate) fn reference(&self) -> &AdditiveFieldGaugeReference {
        &self.reference
    }
}

/// A rejected additive gauge value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum GaugeError {
    /// The selected representative was NaN or infinite.
    NonFiniteFieldValue,
}

impl fmt::Display for GaugeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteFieldValue => formatter.write_str("gauge field value is not finite"),
        }
    }
}

impl Error for GaugeError {}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AdditiveFieldGaugeReference {
    Point(Point3),
    LevelSet(GroupId),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum SharedLevelSetInput {
    General(SharedLevelSet),
    Horizon(Horizon),
}

impl SharedLevelSetInput {
    pub(crate) fn group_id(&self) -> &GroupId {
        match self {
            Self::General(group) => group.group_id(),
            Self::Horizon(group) => group.group_id(),
        }
    }

    pub(crate) fn members(&self) -> &[SharedLevelSetMember] {
        match self {
            Self::General(group) => group.members(),
            Self::Horizon(group) => group.members(),
        }
    }

    pub(crate) fn is_horizon(&self) -> bool {
        matches!(self, Self::Horizon(_))
    }
}

impl private::Sealed for SharedLevelSet {}

impl ProblemInput for SharedLevelSet {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_shared_level_set(SharedLevelSetInput::General(self))
    }
}

impl private::Sealed for Horizon {}

impl ProblemInput for Horizon {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_shared_level_set(SharedLevelSetInput::Horizon(self))
    }
}

impl private::Sealed for AdditiveFieldGauge {}

impl ProblemInput for AdditiveFieldGauge {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_additive_field_gauge(self)
    }
}

impl private::Sealed for FieldValueBound {}

impl ProblemInput for FieldValueBound {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_field_value_bound(self)
    }
}

impl private::Sealed for DirectionalDerivativeInterval {}

impl ProblemInput for DirectionalDerivativeInterval {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_directional_derivative_interval(self)
    }
}

impl private::Sealed for FieldSeparationInterval {}

impl ProblemInput for FieldSeparationInterval {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_field_separation_interval(self)
    }
}

impl private::Sealed for PointToLevelSetRelation {}

impl ProblemInput for PointToLevelSetRelation {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_point_to_level_set_relation(self)
    }
}

impl private::Sealed for YoungerThan {}

impl ProblemInput for YoungerThan {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_shared_level_set_relation(self.0)
    }
}

impl private::Sealed for OlderThan {}

impl ProblemInput for OlderThan {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_shared_level_set_relation(self.0)
    }
}

impl private::Sealed for FieldLevelOrder {}

impl ProblemInput for FieldLevelOrder {
    fn add_to(self, builder: &mut ProblemBuilder) -> Result<(), AddError> {
        builder.add_shared_level_set_relation(self.0)
    }
}
