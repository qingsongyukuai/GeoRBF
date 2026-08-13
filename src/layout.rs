//! Deterministic row/column layouts for the five frozen Surfe models.
//!
//! This module describes indices only. Numerical kernel, polynomial, right-hand
//! side, and smoothing assembly remains outside the T16 boundary.
//!
//! Sources:
//! - `surfe_lib/single_surface.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/lajaunie.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/stratigraphic_surfaces.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/continuous_property.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/vector_field.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`

use std::fmt::{self, Write};

use crate::{model, Axis, Constraints, Error, InternalParameters, ModelType, Parameters};

/// Half-open row/column interval.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct IndexRange {
    start: usize,
    end: usize,
}

impl IndexRange {
    pub const fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub const fn start(self) -> usize {
        self.start
    }

    pub const fn end(self) -> usize {
        self.end
    }

    pub const fn len(self) -> usize {
        self.end - self.start
    }

    pub const fn is_empty(self) -> bool {
        self.start == self.end
    }

    pub const fn contains(self, index: usize) -> bool {
        self.start <= index && index < self.end
    }
}

impl fmt::Display for IndexRange {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}..{}", self.start, self.end)
    }
}

/// Raw category counts supplied to a layout, including categories a model
/// deliberately ignores.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct SourceConstraintCounts {
    pub inequalities: usize,
    pub interfaces: usize,
    pub planars: usize,
    pub tangents: usize,
}

impl SourceConstraintCounts {
    pub(crate) fn from_constraints(constraints: &Constraints) -> Self {
        Self {
            inequalities: constraints.inequalities.len(),
            interfaces: constraints.interfaces.len(),
            planars: constraints.planars.len(),
            tangents: constraints.tangents.len(),
        }
    }
}

/// Stable reference to an input point used by a difference degree of freedom.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum LayoutPointRef {
    Interface(usize),
    Inequality(usize),
}

/// The three source-level meanings of Surfe difference rows, with the two
/// inequality orientations kept distinct.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DifferenceKind {
    SameLevelInterface,
    SequencedInterfaces,
    InequalityBelowUpperInterface,
    InequalityAboveLowerInterface,
}

/// One matrix row/column label.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum LayoutDof {
    InequalityValue {
        index: usize,
    },
    InterfaceValue {
        index: usize,
    },
    Difference {
        kind: DifferenceKind,
        positive: LayoutPointRef,
        negative: LayoutPointRef,
    },
    PlanarDerivative {
        index: usize,
        axis: Axis,
    },
    Tangent {
        index: usize,
    },
    PolynomialTerm {
        index: usize,
    },
}

/// Solver meaning of a row/column label.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum LayoutRole {
    Equality,
    Inequality,
    Bounded,
    Polynomial,
}

impl LayoutRole {
    const fn snapshot_name(self) -> &'static str {
        match self {
            Self::Equality => "equality",
            Self::Inequality => "inequality",
            Self::Bounded => "bounded",
            Self::Polynomial => "polynomial",
        }
    }
}

/// Named contiguous section in one model's fixed row/column order.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum LayoutSectionKind {
    InequalityValues,
    InterfaceValues,
    SequencedInterfaceDifferences,
    SequencedInequalityDifferences,
    SameLevelDifferences,
    PlanarDerivatives,
    Tangents,
    Polynomial,
}

impl LayoutSectionKind {
    const fn snapshot_name(self) -> &'static str {
        match self {
            Self::InequalityValues => "inequality_values",
            Self::InterfaceValues => "interface_values",
            Self::SequencedInterfaceDifferences => "sequenced_interface_differences",
            Self::SequencedInequalityDifferences => "sequenced_inequality_differences",
            Self::SameLevelDifferences => "same_level_differences",
            Self::PlanarDerivatives => "planar_derivatives",
            Self::Tangents => "tangents",
            Self::Polynomial => "polynomial",
        }
    }
}

/// One named contiguous matrix section.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct LayoutSection {
    kind: LayoutSectionKind,
    range: IndexRange,
}

impl LayoutSection {
    pub const fn kind(self) -> LayoutSectionKind {
        self.kind
    }

    pub const fn range(self) -> IndexRange {
        self.range
    }
}

/// Equality, inequality, bounded, and polynomial portions of a layout.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct LayoutPartitions {
    equality: IndexRange,
    inequality: IndexRange,
    bounded: IndexRange,
    polynomial: IndexRange,
}

impl LayoutPartitions {
    pub(crate) const fn new(
        equality: IndexRange,
        inequality: IndexRange,
        bounded: IndexRange,
        polynomial: IndexRange,
    ) -> Self {
        Self {
            equality,
            inequality,
            bounded,
            polynomial,
        }
    }

    pub const fn equality(self) -> IndexRange {
        self.equality
    }

    pub const fn inequality(self) -> IndexRange {
        self.inequality
    }

    pub const fn bounded(self) -> IndexRange {
        self.bounded
    }

    pub const fn polynomial(self) -> IndexRange {
        self.polynomial
    }
}

/// Complete deterministic row/column description for one frozen model.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConstraintLayout {
    model: ModelType,
    source_counts: SourceConstraintCounts,
    internal_parameters: InternalParameters,
    dofs: Vec<LayoutDof>,
    roles: Vec<LayoutRole>,
    sections: Vec<LayoutSection>,
    partitions: LayoutPartitions,
}

impl ConstraintLayout {
    pub(crate) fn new(
        model: ModelType,
        source_counts: SourceConstraintCounts,
        internal_parameters: InternalParameters,
        dofs: Vec<LayoutDof>,
        roles: Vec<LayoutRole>,
        sections: Vec<LayoutSection>,
        partitions: LayoutPartitions,
    ) -> Self {
        debug_assert_eq!(dofs.len(), roles.len());
        Self {
            model,
            source_counts,
            internal_parameters,
            dofs,
            roles,
            sections,
            partitions,
        }
    }

    pub const fn model(&self) -> ModelType {
        self.model
    }

    pub const fn source_counts(&self) -> SourceConstraintCounts {
        self.source_counts
    }

    pub const fn internal_parameters(&self) -> &InternalParameters {
        &self.internal_parameters
    }

    pub fn dofs(&self) -> &[LayoutDof] {
        &self.dofs
    }

    pub fn dof(&self, index: usize) -> Option<&LayoutDof> {
        self.dofs.get(index)
    }

    pub fn role(&self, index: usize) -> Option<LayoutRole> {
        self.roles.get(index).copied()
    }

    pub fn index_of(&self, dof: &LayoutDof) -> Option<usize> {
        self.dofs.iter().position(|candidate| candidate == dof)
    }

    pub fn sections(&self) -> &[LayoutSection] {
        &self.sections
    }

    pub fn section(&self, kind: LayoutSectionKind) -> Option<IndexRange> {
        self.sections
            .iter()
            .find_map(|section| (section.kind == kind).then_some(section.range))
    }

    pub const fn partitions(&self) -> LayoutPartitions {
        self.partitions
    }

    pub fn matrix_size(&self) -> usize {
        self.dofs.len()
    }

    pub fn polynomial_dof_count(&self) -> usize {
        self.partitions.polynomial.len()
    }

    pub fn constraint_dof_count(&self) -> usize {
        self.matrix_size() - self.polynomial_dof_count()
    }

    /// Exact, stable, human-readable evidence for reviews and fixtures.
    pub fn debug_snapshot(&self) -> String {
        let mut snapshot = String::new();
        let internal = &self.internal_parameters;
        writeln!(
            snapshot,
            "model={} solver={:?} modified={} restricted={}",
            self.model, internal.problem_type, internal.modified_basis, internal.restricted_range
        )
        .expect("writing to String cannot fail");
        writeln!(
            snapshot,
            "source inequality={} interface={} planar={} tangent={}",
            self.source_counts.inequalities,
            self.source_counts.interfaces,
            self.source_counts.planars,
            self.source_counts.tangents
        )
        .expect("writing to String cannot fail");
        writeln!(
            snapshot,
            "internal inequality={} interface={} planar={} tangent={} constraints={} equality={} polynomial={}",
            internal.n_inequality,
            internal.n_interface,
            internal.n_planar,
            internal.n_tangent,
            internal.n_constraints,
            internal.n_equality,
            internal.n_poly_terms
        )
        .expect("writing to String cannot fail");
        writeln!(
            snapshot,
            "matrix size={} equality={} inequality={} bounded={} polynomial={}",
            self.matrix_size(),
            self.partitions.equality,
            self.partitions.inequality,
            self.partitions.bounded,
            self.partitions.polynomial
        )
        .expect("writing to String cannot fail");
        for section in &self.sections {
            writeln!(
                snapshot,
                "section {}={}",
                section.kind.snapshot_name(),
                section.range
            )
            .expect("writing to String cannot fail");
        }
        for (index, (dof, role)) in self.dofs.iter().zip(&self.roles).enumerate() {
            write!(
                snapshot,
                "{} {} {}",
                index,
                role.snapshot_name(),
                DofSnapshot(dof)
            )
            .expect("writing to String cannot fail");
            if index + 1 != self.dofs.len() {
                snapshot.push('\n');
            }
        }
        snapshot
    }
}

struct DofSnapshot<'a>(&'a LayoutDof);

impl fmt::Display for DofSnapshot<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            LayoutDof::InequalityValue { index } => write!(formatter, "inequality[{index}]"),
            LayoutDof::InterfaceValue { index } => write!(formatter, "interface[{index}]"),
            LayoutDof::Difference {
                kind,
                positive,
                negative,
            } => write!(
                formatter,
                "difference.{}({}-{})",
                difference_name(*kind),
                PointRefSnapshot(*positive),
                PointRefSnapshot(*negative)
            ),
            LayoutDof::PlanarDerivative { index, axis } => {
                write!(formatter, "planar[{index}].d{}", axis_name(*axis))
            }
            LayoutDof::Tangent { index } => write!(formatter, "tangent[{index}]"),
            LayoutDof::PolynomialTerm { index } => write!(formatter, "polynomial[{index}]"),
        }
    }
}

struct PointRefSnapshot(LayoutPointRef);

impl fmt::Display for PointRefSnapshot {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            LayoutPointRef::Interface(index) => write!(formatter, "interface[{index}]"),
            LayoutPointRef::Inequality(index) => write!(formatter, "inequality[{index}]"),
        }
    }
}

const fn difference_name(kind: DifferenceKind) -> &'static str {
    match kind {
        DifferenceKind::SameLevelInterface => "same_level",
        DifferenceKind::SequencedInterfaces => "sequenced_interfaces",
        DifferenceKind::InequalityBelowUpperInterface => "upper_minus_inequality",
        DifferenceKind::InequalityAboveLowerInterface => "inequality_minus_lower",
    }
}

const fn axis_name(axis: Axis) -> &'static str {
    match axis {
        Axis::X => "x",
        Axis::Y => "y",
        Axis::Z => "z",
    }
}

/// Build the fixed row/column layout for exactly one model.
///
/// The supplied constraints are consumed in their current order. The public
/// Surfe pipeline removes collocations first; callers reproducing that pipeline
/// should call [`Constraints::remove_collocated`] before this function.
pub fn constraint_layout(
    model: ModelType,
    constraints: &Constraints,
    parameters: &Parameters,
) -> Result<ConstraintLayout, Error> {
    match model {
        ModelType::SingleSurface => model::single_surface::layout::build(constraints, parameters),
        ModelType::LajaunieApproach => model::lajaunie::layout::build(constraints, parameters),
        ModelType::StratigraphicHorizons => {
            model::stratigraphic::layout::build(constraints, parameters)
        }
        ModelType::ContinuousProperty => model::continuous_property::layout::build(constraints),
        ModelType::VectorField => model::vector_field::layout::build(constraints),
    }
}

pub(crate) fn section(kind: LayoutSectionKind, start: usize, end: usize) -> LayoutSection {
    LayoutSection {
        kind,
        range: IndexRange::new(start, end),
    }
}

pub(crate) fn declared_polynomial_term_count(order: i32) -> Result<usize, Error> {
    let m = i64::from(order) + 1;
    let terms = m
        .checked_mul(m + 1)
        .and_then(|value| value.checked_mul(m + 2))
        .map(|value| value / 6)
        .ok_or(Error::InterpolationMatrixFailure)?;
    usize::try_from(terms).map_err(|_| Error::InterpolationMatrixFailure)
}

pub(crate) fn planar_dofs(count: usize) -> Vec<LayoutDof> {
    let mut dofs = Vec::with_capacity(count.saturating_mul(3));
    for index in 0..count {
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            dofs.push(LayoutDof::PlanarDerivative { index, axis });
        }
    }
    dofs
}

pub(crate) fn tangent_dofs(count: usize) -> Vec<LayoutDof> {
    (0..count)
        .map(|index| LayoutDof::Tangent { index })
        .collect()
}

pub(crate) fn polynomial_dofs(count: usize) -> Vec<LayoutDof> {
    (0..count)
        .map(|index| LayoutDof::PolynomialTerm { index })
        .collect()
}

pub(crate) fn append_section(
    all_dofs: &mut Vec<LayoutDof>,
    all_roles: &mut Vec<LayoutRole>,
    sections: &mut Vec<LayoutSection>,
    kind: LayoutSectionKind,
    dofs: impl IntoIterator<Item = LayoutDof>,
    role: LayoutRole,
) {
    let start = all_dofs.len();
    for dof in dofs {
        all_dofs.push(dof);
        all_roles.push(role);
    }
    sections.push(section(kind, start, all_dofs.len()));
}
