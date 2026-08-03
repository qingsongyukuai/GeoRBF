use std::cmp::Ordering;
use std::fmt;

use crate::math::canonical_zero;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum FunctionalDimension {
    FieldValue,
    FieldValuePerLength,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FunctionalTerm {
    support: [f64; 3],
    value_coefficient: f64,
    gradient_coefficient: [f64; 3],
}

impl FunctionalTerm {
    pub(crate) fn new(
        support: [f64; 3],
        value_coefficient: f64,
        gradient_coefficient: [f64; 3],
    ) -> Self {
        Self {
            support: support.map(canonical_zero),
            value_coefficient: canonical_zero(value_coefficient),
            gradient_coefficient: gradient_coefficient.map(canonical_zero),
        }
    }

    pub(crate) fn support(self) -> [f64; 3] {
        self.support
    }

    pub(crate) fn value_coefficient(self) -> f64 {
        self.value_coefficient
    }

    pub(crate) fn gradient_coefficient(self) -> [f64; 3] {
        self.gradient_coefficient
    }

    fn is_zero(self) -> bool {
        self.value_coefficient == 0.0
            && self
                .gradient_coefficient
                .iter()
                .all(|coefficient| *coefficient == 0.0)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalFunctional {
    dimension: FunctionalDimension,
    terms: Vec<FunctionalTerm>,
}

impl CanonicalFunctional {
    pub(crate) fn new(
        dimension: FunctionalDimension,
        mut terms: Vec<FunctionalTerm>,
    ) -> Result<Self, FunctionalError> {
        validate_finite(&terms)?;
        terms.sort_by(compare_terms);

        let mut merged: Vec<FunctionalTerm> = Vec::with_capacity(terms.len());
        for term in terms {
            let has_same_support = merged
                .last()
                .is_some_and(|previous| previous.support == term.support);
            if has_same_support {
                let previous = merged
                    .last_mut()
                    .expect("the support match proved a last term");
                previous.value_coefficient =
                    canonical_zero(previous.value_coefficient + term.value_coefficient);
                for axis in 0..3 {
                    previous.gradient_coefficient[axis] = canonical_zero(
                        previous.gradient_coefficient[axis] + term.gradient_coefficient[axis],
                    );
                }
                if !previous.value_coefficient.is_finite()
                    || previous
                        .gradient_coefficient
                        .iter()
                        .any(|coefficient| !coefficient.is_finite())
                {
                    return Err(FunctionalError::NonFiniteMergedCoefficient {
                        support: previous.support,
                    });
                }
                continue;
            }
            merged.push(term);
        }
        merged.retain(|term| !term.is_zero());
        if merged.is_empty() {
            return Err(FunctionalError::ZeroFunctional);
        }

        Ok(Self {
            dimension,
            terms: merged,
        })
    }

    pub(crate) fn dimension(&self) -> FunctionalDimension {
        self.dimension
    }

    pub(crate) fn terms(&self) -> &[FunctionalTerm] {
        &self.terms
    }

    pub(crate) fn evaluate_affine(&self, constant: f64, gradient: [f64; 3]) -> f64 {
        self.terms
            .iter()
            .map(|term| {
                term.value_coefficient
                    * (constant
                        + term.support[0] * gradient[0]
                        + term.support[1] * gradient[1]
                        + term.support[2] * gradient[2])
                    + term.gradient_coefficient[0] * gradient[0]
                    + term.gradient_coefficient[1] * gradient[1]
                    + term.gradient_coefficient[2] * gradient[2]
            })
            .sum()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum FunctionalError {
    NonFiniteSupport { term: usize, axis: usize },
    NonFiniteValueCoefficient { term: usize },
    NonFiniteGradientCoefficient { term: usize, axis: usize },
    NonFiniteMergedCoefficient { support: [f64; 3] },
    ZeroFunctional,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceId(Box<str>);

impl SourceId {
    /// Owns a caller-supplied stable source identity.
    pub fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    /// Returns the caller-supplied identity.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for SourceId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GroupId(Box<str>);

impl GroupId {
    /// Owns a caller-supplied stable group identity.
    pub fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    /// Returns the caller-supplied identity.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for GroupId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct RelationId(Box<str>);

impl RelationId {
    pub(crate) fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct ResidualId(Box<str>);

impl ResidualId {
    pub(crate) fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct DerivedBlockId(Box<str>);

impl DerivedBlockId {
    pub(crate) fn from_residual(residual: &ResidualId) -> Self {
        Self(format!("{}/derived/equality-block", residual.as_str()).into())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct DerivedRowId(Box<str>);

impl DerivedRowId {
    pub(crate) fn from_residual(residual: &ResidualId) -> Self {
        Self(format!("{}/derived/equality-row", residual.as_str()).into())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct DerivedColumnId(Box<str>);

impl DerivedColumnId {
    pub(crate) fn from_residual(residual: &ResidualId) -> Self {
        Self(format!("{}/derived/representer-column", residual.as_str()).into())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SemanticRolePath(Box<str>);

impl SemanticRolePath {
    pub(crate) fn new(value: impl Into<Box<str>>) -> Self {
        Self(value.into())
    }

    /// Returns the stable semantic component path.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for SemanticRolePath {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct UsageProvenance {
    source: SourceId,
    group: Option<GroupId>,
    relation: RelationId,
    residual: ResidualId,
    semantic_role: SemanticRolePath,
}

impl UsageProvenance {
    pub(crate) fn new(
        source: SourceId,
        group: Option<GroupId>,
        relation: RelationId,
        residual: ResidualId,
        semantic_role: SemanticRolePath,
    ) -> Self {
        Self {
            source,
            group,
            relation,
            residual,
            semantic_role,
        }
    }

    pub(crate) fn source(&self) -> &SourceId {
        &self.source
    }

    pub(crate) fn group(&self) -> Option<&GroupId> {
        self.group.as_ref()
    }

    pub(crate) fn relation(&self) -> &RelationId {
        &self.relation
    }

    pub(crate) fn residual(&self) -> &ResidualId {
        &self.residual
    }

    pub(crate) fn semantic_role(&self) -> &SemanticRolePath {
        &self.semantic_role
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FunctionalUse {
    functional: CanonicalFunctional,
    provenance: UsageProvenance,
}

impl FunctionalUse {
    pub(crate) fn new(functional: CanonicalFunctional, provenance: UsageProvenance) -> Self {
        Self {
            functional,
            provenance,
        }
    }

    pub(crate) fn functional(&self) -> &CanonicalFunctional {
        &self.functional
    }

    pub(crate) fn provenance(&self) -> &UsageProvenance {
        &self.provenance
    }
}

fn validate_finite(terms: &[FunctionalTerm]) -> Result<(), FunctionalError> {
    for (term_index, term) in terms.iter().enumerate() {
        for (axis, coordinate) in term.support.iter().enumerate() {
            if !coordinate.is_finite() {
                return Err(FunctionalError::NonFiniteSupport {
                    term: term_index,
                    axis,
                });
            }
        }
        if !term.value_coefficient.is_finite() {
            return Err(FunctionalError::NonFiniteValueCoefficient { term: term_index });
        }
        for (axis, coefficient) in term.gradient_coefficient.iter().enumerate() {
            if !coefficient.is_finite() {
                return Err(FunctionalError::NonFiniteGradientCoefficient {
                    term: term_index,
                    axis,
                });
            }
        }
    }
    Ok(())
}

fn compare_terms(left: &FunctionalTerm, right: &FunctionalTerm) -> Ordering {
    left.support
        .iter()
        .chain(std::iter::once(&left.value_coefficient))
        .chain(left.gradient_coefficient.iter())
        .zip(
            right
                .support
                .iter()
                .chain(std::iter::once(&right.value_coefficient))
                .chain(right.gradient_coefficient.iter()),
        )
        .find_map(|(left, right)| {
            let ordering = left.total_cmp(right);
            (ordering != Ordering::Equal).then_some(ordering)
        })
        .unwrap_or(Ordering::Equal)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn term(
        support: [f64; 3],
        value_coefficient: f64,
        gradient_coefficient: [f64; 3],
    ) -> FunctionalTerm {
        FunctionalTerm::new(support, value_coefficient, gradient_coefficient)
    }

    #[test]
    fn canonical_functional_merges_and_sorts_exact_supports() {
        let functional = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![
                term([2.0, 0.0, 1.0], 1.0, [0.0, 2.0, -0.0]),
                term([-1.0, 3.0, 0.0], 4.0, [1.0, 0.0, 0.0]),
                term([2.0, -0.0, 1.0], -1.0, [-0.0, -2.0, 0.0]),
                term([0.0, 0.0, 0.0], -0.0, [0.0; 3]),
                term([-1.0, 3.0, -0.0], -1.5, [0.0, 0.5, 0.0]),
            ],
        )
        .expect("a nonzero functional should canonicalize");

        assert_eq!(
            functional.terms(),
            &[term([-1.0, 3.0, 0.0], 2.5, [1.0, 0.5, 0.0])]
        );
        assert!(
            functional
                .terms()
                .iter()
                .flat_map(|term| {
                    std::iter::once(term.value_coefficient()).chain(term.gradient_coefficient())
                })
                .filter(|value| *value == 0.0)
                .all(|value| !value.is_sign_negative())
        );
    }

    #[test]
    fn canonical_functional_is_permutation_invariant_but_does_not_merge_nearby_supports() {
        let one = term([0.0, 0.0, 0.0], 1.0, [0.0; 3]);
        let nearby = term([f64::EPSILON, 0.0, 0.0], -1.0, [0.0; 3]);
        let forward = CanonicalFunctional::new(FunctionalDimension::FieldValue, vec![one, nearby])
            .expect("distinct supports keep the functional nonzero");
        let reversed = CanonicalFunctional::new(FunctionalDimension::FieldValue, vec![nearby, one])
            .expect("input order must not matter");

        assert_eq!(forward, reversed);
        assert_eq!(forward.terms().len(), 2);
    }

    #[test]
    fn zero_functional_is_rejected_structurally() {
        let failure = CanonicalFunctional::new(
            FunctionalDimension::FieldValuePerLength,
            vec![term([1.0, 2.0, 3.0], 0.0, [-0.0, 0.0, 0.0])],
        )
        .expect_err("zero functionals have no representer");

        assert_eq!(failure, FunctionalError::ZeroFunctional);
    }

    #[test]
    fn dimension_belongs_to_the_functional_and_provenance_belongs_to_its_usage_edge() {
        let functional = CanonicalFunctional::new(
            FunctionalDimension::FieldValue,
            vec![term([1.0, 2.0, 3.0], 2.0, [1.0, -0.5, 0.25])],
        )
        .expect("fixture functional is valid");
        let first = FunctionalUse::new(
            functional.clone(),
            UsageProvenance::new(
                SourceId::new("source-a"),
                Some(GroupId::new("group-a")),
                RelationId::new("relation-a"),
                ResidualId::new("residual-a"),
                SemanticRolePath::new("hard-equality/left"),
            ),
        );
        let second = FunctionalUse::new(
            functional.clone(),
            UsageProvenance::new(
                SourceId::new("source-b"),
                Some(GroupId::new("group-b")),
                RelationId::new("relation-b"),
                ResidualId::new("residual-b"),
                SemanticRolePath::new("hard-equality/right"),
            ),
        );

        assert_eq!(first.functional(), second.functional());
        assert_ne!(first.provenance(), second.provenance());
        assert_eq!(functional.dimension(), FunctionalDimension::FieldValue);
        assert_eq!(first.provenance().source(), &SourceId::new("source-a"));
        assert_eq!(first.provenance().group(), Some(&GroupId::new("group-a")));
        assert_eq!(
            first.provenance().relation(),
            &RelationId::new("relation-a")
        );
        assert_eq!(
            first.provenance().residual(),
            &ResidualId::new("residual-a")
        );
        assert_eq!(
            first.provenance().semantic_role(),
            &SemanticRolePath::new("hard-equality/left")
        );
    }

    #[test]
    fn four_term_difference_and_normal_tangent_functionals_have_analytic_truth() {
        use crate::cubic::{CubicKernel, GlobalAnisotropyMetric};

        let difference = |left: f64, right: f64| {
            CanonicalFunctional::new(
                FunctionalDimension::FieldValue,
                vec![
                    term([left, 0.0, 0.0], 1.0, [0.0; 3]),
                    term([right, 0.0, 0.0], -1.0, [0.0; 3]),
                ],
            )
            .expect("a difference with distinct supports is nonzero")
        };
        assert_eq!(
            CubicKernel::pairing(
                &difference(0.0, 1.0),
                &difference(3.0, 5.0),
                &GlobalAnisotropyMetric::identity(),
            ),
            -42.0
        );

        let normal = CanonicalFunctional::new(
            FunctionalDimension::FieldValuePerLength,
            vec![term([0.0; 3], 0.0, [1.0, 0.0, 0.0])],
        )
        .expect("normal projection functional is nonzero");
        let tangent = CanonicalFunctional::new(
            FunctionalDimension::FieldValuePerLength,
            vec![term([0.0; 3], 0.0, [0.0, 1.0, 0.0])],
        )
        .expect("tangent functional is nonzero");
        assert_eq!(normal.evaluate_affine(1.0, [2.0, 0.0, -0.5]), 2.0);
        assert_eq!(tangent.evaluate_affine(1.0, [2.0, 0.0, -0.5]), 0.0);
    }
}
