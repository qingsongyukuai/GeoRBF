//! Shared level sets, geological horizons, and explicit additive field gauges.

use std::collections::BTreeSet;
use std::error::Error;
use std::fmt;

use crate::functional::{GroupId, SourceId};
use crate::geometry::Point3;
use crate::problem::{AddError, ProblemBuilder, ProblemInput, private};

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
