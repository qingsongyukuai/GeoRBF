//! Solver-independent assembly of canonical Cubic objectives and relations.
//!
//! The representation owns numerical coordinates. This module owns canonical
//! hard/soft semantics and records every source-bearing row once; backend
//! adapters only realize one of the coordinate layouts retained by the form.

use std::collections::{BTreeMap, BTreeSet};

use crate::cubic_equality::{
    CanonicalEqualityParticipation, CanonicalHardSourceRecovery, CanonicalInequalitySense,
    CanonicalRelationToleranceEvidence, CanonicalSoftLoss, CanonicalSoftResidualBlockKind,
    CanonicalViolationLoss, CpdEvidence, CubicCanonicalProblem, CubicFunctionalResponse,
    CubicRepresentation, RecoveredHardEquality, RecoveredSoftEquality, RecoveredSoftObjective,
    RepresentationFailure, SemanticLatentCoefficient,
};
use crate::functional::{
    CanonicalFunctional, DerivedBlockId, DerivedColumnId, DerivedRowId, FunctionalDimension,
    GroupId, ResidualId, SourceId, UsageProvenance,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AllSourceRecoveryLedger {
    pub(crate) canonical_hard_relations: usize,
    pub(crate) canonical_soft_relations: usize,
    pub(crate) participating_sources: Vec<SourceId>,
    pub(crate) recovered_sources: Vec<SourceId>,
    pub(crate) representers: usize,
    pub(crate) solver_relation_rows: usize,
    pub(crate) recovery_edges: usize,
    pub(crate) verified: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecoveredAffineRelation {
    pub(crate) canonical_index: usize,
    pub(crate) provenances: Vec<UsageProvenance>,
    pub(crate) violation_loss: Option<CanonicalViolationLoss>,
    pub(crate) violation: f64,
    pub(crate) tolerance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CubicFieldCoordinateLayout {
    Standard,
    Quotient,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CubicSolverVariableLayout {
    pub(crate) field: usize,
    pub(crate) polynomial: usize,
    pub(crate) semantic_latents: usize,
    pub(crate) auxiliary: usize,
}

impl CubicSolverVariableLayout {
    pub(crate) fn variables(self) -> Option<usize> {
        self.field
            .checked_add(self.polynomial)?
            .checked_add(self.semantic_latents)?
            .checked_add(self.auxiliary)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSolverRow {
    pub(crate) canonical_index: usize,
    pub(crate) response: Option<CubicFunctionalResponse>,
    pub(crate) latent_coefficients: Vec<SemanticLatentCoefficient>,
    pub(crate) provenance: UsageProvenance,
    pub(crate) source_provenances: Vec<UsageProvenance>,
    pub(crate) derived_block: DerivedBlockId,
    pub(crate) residual: ResidualId,
    pub(crate) derived_row: DerivedRowId,
    pub(crate) derived_column: Option<DerivedColumnId>,
    pub(crate) dimension: FunctionalDimension,
    pub(crate) target: f64,
}

struct CanonicalSolverRowInput<'a> {
    canonical_index: usize,
    functional: Option<&'a CanonicalFunctional>,
    latent_coefficients: &'a [SemanticLatentCoefficient],
    provenance: UsageProvenance,
    source_provenances: Vec<UsageProvenance>,
    dimension: FunctionalDimension,
    target: f64,
}

impl CanonicalSolverRow {
    fn from_canonical(
        representation: &CubicRepresentation,
        input: CanonicalSolverRowInput<'_>,
    ) -> Result<Self, RepresentationFailure> {
        let residual = input.provenance.residual().clone();
        Ok(Self {
            canonical_index: input.canonical_index,
            response: input
                .functional
                .map(|functional| representation.response(functional))
                .transpose()?,
            latent_coefficients: input.latent_coefficients.to_vec(),
            source_provenances: input.source_provenances,
            derived_block: DerivedBlockId::from_residual(&residual),
            residual: residual.clone(),
            derived_row: DerivedRowId::from_residual(&residual),
            derived_column: input
                .functional
                .map(|_| DerivedColumnId::from_residual(&residual)),
            dimension: input.dimension,
            target: input.target,
            provenance: input.provenance,
        })
    }

    pub(crate) fn coefficients(
        &self,
        coordinate_layout: CubicFieldCoordinateLayout,
        variable_layout: CubicSolverVariableLayout,
    ) -> Vec<f64> {
        let variables = variable_layout
            .variables()
            .expect("a validated Cubic solver form has a finite variable count");
        let mut coefficients = vec![0.0; variables];
        if let Some(response) = &self.response {
            let field = match coordinate_layout {
                CubicFieldCoordinateLayout::Standard => &response.standard_field,
                CubicFieldCoordinateLayout::Quotient => &response.quotient_field,
            };
            debug_assert_eq!(field.len(), variable_layout.field);
            coefficients[..variable_layout.field].copy_from_slice(field);
            coefficients[variable_layout.field..variable_layout.field + variable_layout.polynomial]
                .copy_from_slice(&response.polynomial);
        }
        let latent_offset = variable_layout.field + variable_layout.polynomial;
        for term in &self.latent_coefficients {
            coefficients[latent_offset + term.latent] = term.coefficient;
        }
        coefficients
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardSolverRow {
    pub(crate) row: CanonicalSolverRow,
    pub(crate) participation: CanonicalEqualityParticipation,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardRecoveryRelation {
    pub(crate) canonical_index: usize,
    pub(crate) provenance: UsageProvenance,
    pub(crate) target: f64,
    pub(crate) relation_to_canonical: f64,
    /// `(retained solver row, coefficient)` entries reconstructing this row.
    pub(crate) coefficients: Vec<(usize, f64)>,
    pub(crate) complete_affine_verified: bool,
}

#[derive(Debug, Clone, PartialEq)]
struct CanonicalHardRowRecovery {
    canonical_index: usize,
    coefficients: Vec<(usize, f64)>,
    complete_affine_verified: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardRecoveryGraph {
    pub(crate) retained_rows: Vec<usize>,
    rows: Vec<CanonicalHardRowRecovery>,
    pub(crate) relations: Vec<CanonicalHardRecoveryRelation>,
    pub(crate) conflict_witness: Option<CanonicalHardConflictWitness>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardConflictRelation {
    pub(crate) provenance: UsageProvenance,
    pub(crate) multiplier: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalHardConflictWitness {
    pub(crate) relations: Vec<CanonicalHardConflictRelation>,
    pub(crate) sources: Vec<UsageProvenance>,
    pub(crate) canonical_residual: f64,
    pub(crate) separation_margin: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct SymbolicAffineRow {
    coefficients: BTreeMap<usize, f64>,
    target: f64,
    source_recoveries: Vec<CanonicalHardSourceRecovery>,
    solver_constraint: bool,
}

impl SymbolicAffineRow {
    #[cfg(test)]
    fn new(coefficients: Vec<f64>, target: f64) -> Self {
        Self {
            coefficients: coefficients
                .into_iter()
                .enumerate()
                .filter(|(_, coefficient)| *coefficient != 0.0)
                .collect(),
            target,
            source_recoveries: Vec::new(),
            solver_constraint: true,
        }
    }
}

impl CanonicalHardRecoveryGraph {
    fn verifies(&self, rows: &[SymbolicAffineRow]) -> bool {
        self.rows.len() == rows.len()
            && self
                .rows
                .iter()
                .enumerate()
                .all(|(canonical_index, recovery)| {
                    let row = &rows[canonical_index];
                    recovery.canonical_index == canonical_index
                        && if recovery.complete_affine_verified {
                            exactly_reconstructs(
                                row,
                                &recovery.coefficients,
                                &self.retained_rows,
                                rows,
                            )
                        } else {
                            !row.solver_constraint
                                && recovery
                                    .coefficients
                                    .iter()
                                    .all(|(solver_row, coefficient)| {
                                        *solver_row < self.retained_rows.len()
                                            && coefficient.is_finite()
                                    })
                        }
                })
            && self.relations.len()
                == rows
                    .iter()
                    .map(|row| row.source_recoveries.len())
                    .sum::<usize>()
            && self.relations.iter().all(|relation| {
                self.rows.get(relation.canonical_index).is_some_and(|row| {
                    rows[relation.canonical_index]
                        .source_recoveries
                        .iter()
                        .find(|source| {
                            source.provenance == relation.provenance
                                && source.target == relation.target
                        })
                        .is_some_and(|source| {
                            relation.coefficients
                                == row
                                    .coefficients
                                    .iter()
                                    .map(|(solver_row, coefficient)| {
                                        (*solver_row, source.coefficient * coefficient)
                                    })
                                    .collect::<Vec<_>>()
                                && relation.complete_affine_verified
                                    == (row.complete_affine_verified
                                        && exactly_reconstructs_source(
                                            &rows[relation.canonical_index],
                                            source,
                                            &relation.coefficients,
                                            &self.retained_rows,
                                            rows,
                                        ))
                                && relation.relation_to_canonical == source.coefficient
                        })
                })
            })
    }

    fn verifies_conflict_witness(&self, rows: &[SymbolicAffineRow]) -> bool {
        let Some(witness) = &self.conflict_witness else {
            return true;
        };
        if witness.relations.is_empty()
            || witness.canonical_residual != 0.0
            || !witness.separation_margin.is_finite()
            || witness.separation_margin <= 0.0
        {
            return false;
        }
        let recovered = witness
            .relations
            .iter()
            .filter_map(|relation| {
                rows.iter().find_map(|row| {
                    row.source_recoveries
                        .iter()
                        .find(|source| source.provenance == relation.provenance)
                        .map(|source| (row, source, relation.multiplier))
                })
            })
            .collect::<Vec<_>>();
        if recovered.len() != witness.relations.len() {
            return false;
        }
        let mut coefficients = BTreeMap::<usize, ExactBinaryProductSum>::new();
        for (row, source, multiplier) in &recovered {
            for (column, coefficient) in &row.coefficients {
                let scaled_source =
                    ExactBinaryProductSum::product(source.coefficient, *coefficient);
                let product = scaled_source.multiply_by_f64(*multiplier);
                let sum = coefficients
                    .remove(column)
                    .unwrap_or_else(ExactBinaryProductSum::zero)
                    .add(&product);
                if !sum.magnitude.is_empty() {
                    coefficients.insert(*column, sum);
                }
            }
        }
        let recomputed_margin = conflict_separation_margin(
            recovered
                .iter()
                .map(|(_, source, multiplier)| (*multiplier, source.target)),
        );
        let mut recovered_sources = recovered
            .iter()
            .flat_map(|(row, _, _)| {
                row.source_recoveries
                    .iter()
                    .map(|source| source.provenance.clone())
            })
            .collect::<Vec<_>>();
        normalize_conflict_sources(&mut recovered_sources);
        coefficients.is_empty()
            && !exact_linear_combination_equals(
                recovered
                    .iter()
                    .map(|(_, source, multiplier)| (*multiplier, source.target)),
                0.0,
            )
            && recomputed_margin
                .is_some_and(|margin| margin.to_bits() == witness.separation_margin.to_bits())
            && recovered_sources == witness.sources
    }
}

#[derive(Debug)]
struct EliminationBasisRow {
    reduced: BTreeMap<usize, f64>,
    retained_expression: BTreeMap<usize, f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ExactBinaryProductSum {
    negative: bool,
    magnitude: Vec<u64>,
}

impl ExactBinaryProductSum {
    fn zero() -> Self {
        Self {
            negative: false,
            magnitude: Vec::new(),
        }
    }

    fn product(left: f64, right: f64) -> Self {
        debug_assert!(left.is_finite() && right.is_finite());
        let (left_negative, left_magnitude) = exact_binary64_magnitude(left);
        let (right_negative, right_magnitude) = exact_binary64_magnitude(right);
        Self::normalized(
            left_negative != right_negative,
            multiply_magnitudes(&left_magnitude, &right_magnitude),
        )
    }

    fn multiply_by_f64(&self, factor: f64) -> Self {
        let (factor_negative, factor_magnitude) = exact_binary64_magnitude(factor);
        Self::normalized(
            self.negative != factor_negative,
            multiply_magnitudes(&self.magnitude, &factor_magnitude),
        )
    }

    fn add(&self, other: &Self) -> Self {
        if self.negative == other.negative {
            return Self::normalized(
                self.negative,
                add_magnitudes(&self.magnitude, &other.magnitude),
            );
        }
        match compare_magnitudes(&self.magnitude, &other.magnitude) {
            std::cmp::Ordering::Greater => Self::normalized(
                self.negative,
                subtract_magnitudes(&self.magnitude, &other.magnitude),
            ),
            std::cmp::Ordering::Less => Self::normalized(
                other.negative,
                subtract_magnitudes(&other.magnitude, &self.magnitude),
            ),
            std::cmp::Ordering::Equal => Self::zero(),
        }
    }

    fn normalized(negative: bool, mut magnitude: Vec<u64>) -> Self {
        while magnitude.last() == Some(&0) {
            magnitude.pop();
        }
        Self {
            negative: negative && !magnitude.is_empty(),
            magnitude,
        }
    }

    fn saturating_f64_magnitude(&self) -> f64 {
        let Some((limb_index, highest_limb)) =
            self.magnitude.iter().copied().enumerate().next_back()
        else {
            return 0.0;
        };
        let highest_bit = 63 - highest_limb.leading_zeros() as i32;
        let exponent = limb_index as i32 * 64 + highest_bit - 2148;
        let significand = highest_limb as f64 / 2.0_f64.powi(highest_bit);
        let magnitude = significand * 2.0_f64.powi(exponent);
        if magnitude.is_infinite() {
            f64::MAX
        } else if magnitude == 0.0 {
            f64::from_bits(1)
        } else {
            magnitude
        }
    }
}

fn exact_binary64_magnitude(value: f64) -> (bool, Vec<u64>) {
    if value == 0.0 {
        return (false, Vec::new());
    }
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let raw_exponent = ((bits >> 52) & 0x7ff) as usize;
    let fraction = bits & ((1_u64 << 52) - 1);
    let (significand, shift) = if raw_exponent == 0 {
        (fraction, 0)
    } else {
        ((1_u64 << 52) | fraction, raw_exponent - 1)
    };
    let limb = shift / 64;
    let bit = shift % 64;
    let mut magnitude = vec![0_u64; limb + 2];
    magnitude[limb] = significand << bit;
    if bit != 0 {
        magnitude[limb + 1] = significand >> (64 - bit);
    }
    while magnitude.last() == Some(&0) {
        magnitude.pop();
    }
    (negative, magnitude)
}

fn compare_magnitudes(left: &[u64], right: &[u64]) -> std::cmp::Ordering {
    left.len()
        .cmp(&right.len())
        .then_with(|| left.iter().rev().cmp(right.iter().rev()))
}

fn add_magnitudes(left: &[u64], right: &[u64]) -> Vec<u64> {
    let mut result = Vec::with_capacity(left.len().max(right.len()) + 1);
    let mut carry = 0_u128;
    for index in 0..left.len().max(right.len()) {
        let sum = u128::from(left.get(index).copied().unwrap_or(0))
            + u128::from(right.get(index).copied().unwrap_or(0))
            + carry;
        result.push(sum as u64);
        carry = sum >> 64;
    }
    if carry != 0 {
        result.push(carry as u64);
    }
    result
}

fn subtract_magnitudes(larger: &[u64], smaller: &[u64]) -> Vec<u64> {
    let mut result = Vec::with_capacity(larger.len());
    let mut borrow = false;
    for (index, larger_limb) in larger.iter().copied().enumerate() {
        let (difference, first_borrow) =
            larger_limb.overflowing_sub(smaller.get(index).copied().unwrap_or(0));
        let (difference, second_borrow) = difference.overflowing_sub(u64::from(borrow));
        result.push(difference);
        borrow = first_borrow || second_borrow;
    }
    debug_assert!(!borrow);
    result
}

fn multiply_magnitudes(left: &[u64], right: &[u64]) -> Vec<u64> {
    if left.is_empty() || right.is_empty() {
        return Vec::new();
    }
    let mut result = vec![0_u64; left.len() + right.len()];
    for (left_index, left_limb) in left.iter().copied().enumerate() {
        let mut carry = 0_u128;
        for (right_index, right_limb) in right.iter().copied().enumerate() {
            let index = left_index + right_index;
            let product =
                u128::from(left_limb) * u128::from(right_limb) + u128::from(result[index]) + carry;
            result[index] = product as u64;
            carry = product >> 64;
        }
        let mut index = left_index + right.len();
        while carry != 0 {
            let sum = u128::from(result[index]) + carry;
            result[index] = sum as u64;
            carry = sum >> 64;
            index += 1;
            if index == result.len() && carry != 0 {
                result.push(0);
            }
        }
    }
    result
}

fn exact_linear_combination_equals(terms: impl Iterator<Item = (f64, f64)>, expected: f64) -> bool {
    terms
        .map(|(coefficient, value)| ExactBinaryProductSum::product(coefficient, value))
        .fold(ExactBinaryProductSum::zero(), |sum, product| {
            sum.add(&product)
        })
        == ExactBinaryProductSum::product(1.0, expected)
}

fn exact_sparse_linear_combination(
    coefficients: &[(usize, f64)],
    retained_rows: &[usize],
    rows: &[SymbolicAffineRow],
) -> BTreeMap<usize, ExactBinaryProductSum> {
    let mut reconstructed = BTreeMap::<usize, ExactBinaryProductSum>::new();
    for (solver_row, recovery_coefficient) in coefficients {
        for (column, value) in &rows[retained_rows[*solver_row]].coefficients {
            let product = ExactBinaryProductSum::product(*recovery_coefficient, *value);
            let sum = reconstructed
                .remove(column)
                .unwrap_or_else(ExactBinaryProductSum::zero)
                .add(&product);
            if !sum.magnitude.is_empty() {
                reconstructed.insert(*column, sum);
            }
        }
    }
    reconstructed
}

fn exactly_reconstructs(
    row: &SymbolicAffineRow,
    coefficients: &[(usize, f64)],
    retained_rows: &[usize],
    rows: &[SymbolicAffineRow],
) -> bool {
    exactly_reconstructs_coefficients(row, coefficients, retained_rows, rows)
        && exact_linear_combination_equals(
            coefficients.iter().map(|(solver_row, coefficient)| {
                (*coefficient, rows[retained_rows[*solver_row]].target)
            }),
            row.target,
        )
}

fn exactly_reconstructs_coefficients(
    row: &SymbolicAffineRow,
    coefficients: &[(usize, f64)],
    retained_rows: &[usize],
    rows: &[SymbolicAffineRow],
) -> bool {
    let reconstructed = exact_sparse_linear_combination(coefficients, retained_rows, rows);
    reconstructed.len() == row.coefficients.len()
        && row.coefficients.iter().all(|(column, expected)| {
            reconstructed.get(column) == Some(&ExactBinaryProductSum::product(1.0, *expected))
        })
}

fn exactly_reconstructs_source(
    canonical_row: &SymbolicAffineRow,
    source: &CanonicalHardSourceRecovery,
    coefficients: &[(usize, f64)],
    retained_rows: &[usize],
    rows: &[SymbolicAffineRow],
) -> bool {
    let reconstructed = exact_sparse_linear_combination(coefficients, retained_rows, rows);
    reconstructed.len() == canonical_row.coefficients.len()
        && canonical_row.coefficients.iter().all(|(column, value)| {
            reconstructed.get(column)
                == Some(&ExactBinaryProductSum::product(source.coefficient, *value))
        })
        && exact_linear_combination_equals(
            coefficients.iter().map(|(solver_row, coefficient)| {
                (*coefficient, rows[retained_rows[*solver_row]].target)
            }),
            source.target,
        )
}

fn build_hard_recovery_graph(rows: &[SymbolicAffineRow]) -> CanonicalHardRecoveryGraph {
    let mut retained_rows = Vec::new();
    let mut row_recoveries = Vec::with_capacity(rows.len());
    let mut relations = Vec::new();
    let mut basis = BTreeMap::<usize, EliminationBasisRow>::new();
    let mut conflict_witness = None;

    for (canonical_index, row) in rows.iter().enumerate() {
        let mut reduced = row.coefficients.clone();
        let mut retained_expression = BTreeMap::<usize, f64>::new();
        while let Some((pivot, basis_row)) = reduced
            .keys()
            .find_map(|pivot| basis.get(pivot).map(|row| (*pivot, row)))
        {
            let factor = reduced[&pivot] / basis_row.reduced[&pivot];
            for (column, basis_value) in &basis_row.reduced {
                let value = reduced.get(column).copied().unwrap_or(0.0) - factor * basis_value;
                if value == 0.0 {
                    reduced.remove(column);
                } else {
                    reduced.insert(*column, value);
                }
            }
            for (solver_row, basis_coefficient) in &basis_row.retained_expression {
                let coefficient = retained_expression.get(solver_row).copied().unwrap_or(0.0)
                    + factor * basis_coefficient;
                if coefficient == 0.0 {
                    retained_expression.remove(solver_row);
                } else {
                    retained_expression.insert(*solver_row, coefficient);
                }
            }
        }
        let proposed = retained_expression
            .iter()
            .map(|(row, value)| (*row, *value))
            .collect::<Vec<_>>();
        let left_is_exact_dependency = reduced.is_empty()
            && exactly_reconstructs_coefficients(row, &proposed, &retained_rows, rows);
        let complete_affine_verified =
            left_is_exact_dependency && exactly_reconstructs(row, &proposed, &retained_rows, rows);
        if left_is_exact_dependency && !complete_affine_verified && conflict_witness.is_none() {
            conflict_witness = hard_dependency_conflict_witness(
                canonical_index,
                row,
                &proposed,
                &retained_rows,
                rows,
            );
        }
        let (coefficients, retained) = if complete_affine_verified || !row.solver_constraint {
            (proposed.clone(), false)
        } else {
            let solver_row = retained_rows.len();
            retained_rows.push(canonical_index);
            let mut basis_expression = retained_expression
                .into_iter()
                .map(|(row, coefficient)| (row, -coefficient))
                .collect::<BTreeMap<_, _>>();
            basis_expression.insert(solver_row, 1.0);
            if let Some(pivot) = reduced.keys().next().copied() {
                basis.insert(
                    pivot,
                    EliminationBasisRow {
                        reduced,
                        retained_expression: basis_expression,
                    },
                );
            }
            (vec![(solver_row, 1.0)], true)
        };
        row_recoveries.push(CanonicalHardRowRecovery {
            canonical_index,
            coefficients: coefficients.clone(),
            complete_affine_verified: complete_affine_verified || retained,
        });
        let complete_affine_verified = row_recoveries
            .last()
            .expect("the recovery row was just appended")
            .complete_affine_verified;
        relations.extend(row.source_recoveries.iter().map(|source| {
            let source_coefficients = coefficients
                .iter()
                .map(|(solver_row, coefficient)| (*solver_row, source.coefficient * coefficient))
                .collect::<Vec<_>>();
            CanonicalHardRecoveryRelation {
                canonical_index,
                provenance: source.provenance.clone(),
                target: source.target,
                relation_to_canonical: source.coefficient,
                complete_affine_verified: complete_affine_verified
                    && exactly_reconstructs_source(
                        row,
                        source,
                        &source_coefficients,
                        &retained_rows,
                        rows,
                    ),
                coefficients: source_coefficients,
            }
        }));
    }
    CanonicalHardRecoveryGraph {
        retained_rows,
        rows: row_recoveries,
        relations,
        conflict_witness,
    }
}

fn hard_dependency_conflict_witness(
    canonical_index: usize,
    row: &SymbolicAffineRow,
    proposed: &[(usize, f64)],
    retained_rows: &[usize],
    rows: &[SymbolicAffineRow],
) -> Option<CanonicalHardConflictWitness> {
    if row.source_recoveries.first()?.coefficient == 0.0 {
        return None;
    }
    let mut canonical_terms = proposed
        .iter()
        .map(|(solver_row, coefficient)| (retained_rows[*solver_row], -*coefficient))
        .collect::<Vec<_>>();
    canonical_terms.push((canonical_index, 1.0));
    let multiplier_scale = canonical_terms
        .iter()
        .map(|(_, multiplier)| multiplier.abs())
        .fold(0.0_f64, f64::max);
    if !multiplier_scale.is_finite() || multiplier_scale == 0.0 {
        return None;
    }
    for (_, multiplier) in &mut canonical_terms {
        *multiplier /= multiplier_scale;
    }
    let mut relations = canonical_terms
        .iter()
        .filter_map(|(index, multiplier)| {
            let source = rows[*index].source_recoveries.first()?;
            (source.coefficient != 0.0).then(|| CanonicalHardConflictRelation {
                provenance: source.provenance.clone(),
                multiplier: *multiplier / source.coefficient,
            })
        })
        .collect::<Vec<_>>();
    relations.sort_by(|left, right| {
        left.provenance
            .source()
            .cmp(right.provenance.source())
            .then_with(|| {
                left.provenance
                    .semantic_role()
                    .cmp(right.provenance.semantic_role())
            })
    });
    let separation_margin = conflict_separation_margin(relations.iter().filter_map(|relation| {
        rows.iter().find_map(|row| {
            row.source_recoveries
                .iter()
                .find(|source| source.provenance == relation.provenance)
                .map(|source| (relation.multiplier, source.target))
        })
    }))?;
    let mut sources = canonical_terms
        .iter()
        .flat_map(|(index, _)| {
            rows[*index]
                .source_recoveries
                .iter()
                .map(|source| source.provenance.clone())
        })
        .collect::<Vec<_>>();
    normalize_conflict_sources(&mut sources);
    Some(CanonicalHardConflictWitness {
        relations,
        sources,
        canonical_residual: 0.0,
        separation_margin,
    })
}

fn normalize_conflict_sources(sources: &mut Vec<UsageProvenance>) {
    sources.sort_by(|left, right| {
        left.source()
            .cmp(right.source())
            .then_with(|| left.semantic_role().cmp(right.semantic_role()))
            .then_with(|| left.groups().cmp(right.groups()))
    });
    sources.dedup_by(|left, right| {
        left.source() == right.source()
            && left.semantic_role() == right.semantic_role()
            && left.groups() == right.groups()
    });
}

fn conflict_separation_margin(terms: impl IntoIterator<Item = (f64, f64)>) -> Option<f64> {
    let exact = terms
        .into_iter()
        .map(|(multiplier, target)| ExactBinaryProductSum::product(multiplier, target))
        .fold(ExactBinaryProductSum::zero(), |sum, product| {
            sum.add(&product)
        });
    let margin = exact.saturating_f64_magnitude();
    (margin > 0.0).then_some(margin)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum SymbolicFieldComponent {
    Value,
    Gradient(usize),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum SymbolicAffineColumn {
    Field {
        dimension: FunctionalDimension,
        support: [u64; 3],
        component: SymbolicFieldComponent,
    },
    SemanticLatent {
        dimension: FunctionalDimension,
        latent: usize,
    },
}

fn symbolic_hard_rows(problem: &CubicCanonicalProblem) -> Vec<SymbolicAffineRow> {
    let columns = problem
        .equalities
        .iter()
        .flat_map(symbolic_columns)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .enumerate()
        .map(|(index, column)| (column, index))
        .collect::<BTreeMap<_, _>>();
    problem
        .equalities
        .iter()
        .map(|equality| {
            let mut coefficients = BTreeMap::new();
            if let Some(field) = equality.field() {
                for term in field.functional().terms() {
                    let support = term.support().map(|coordinate| coordinate.to_bits());
                    if term.value_coefficient() != 0.0 {
                        coefficients.insert(
                            columns[&SymbolicAffineColumn::Field {
                                dimension: equality.dimension(),
                                support,
                                component: SymbolicFieldComponent::Value,
                            }],
                            term.value_coefficient(),
                        );
                    }
                    for (axis, coefficient) in term.gradient_coefficient().into_iter().enumerate() {
                        if coefficient != 0.0 {
                            coefficients.insert(
                                columns[&SymbolicAffineColumn::Field {
                                    dimension: equality.dimension(),
                                    support,
                                    component: SymbolicFieldComponent::Gradient(axis),
                                }],
                                coefficient,
                            );
                        }
                    }
                }
            }
            for term in equality.latent_coefficients() {
                if term.coefficient != 0.0 {
                    coefficients.insert(
                        columns[&SymbolicAffineColumn::SemanticLatent {
                            dimension: equality.dimension(),
                            latent: term.latent,
                        }],
                        term.coefficient,
                    );
                }
            }
            SymbolicAffineRow {
                coefficients,
                target: equality.target(),
                source_recoveries: equality.source_recoveries().to_vec(),
                solver_constraint: equality.participation()
                    == CanonicalEqualityParticipation::SolverConstraint,
            }
        })
        .collect()
}

fn symbolic_columns(
    equality: &crate::cubic_equality::CanonicalHardEquality,
) -> Vec<SymbolicAffineColumn> {
    let mut columns = Vec::new();
    if let Some(field) = equality.field() {
        for term in field.functional().terms() {
            let support = term.support().map(|coordinate| coordinate.to_bits());
            columns.push(SymbolicAffineColumn::Field {
                dimension: equality.dimension(),
                support,
                component: SymbolicFieldComponent::Value,
            });
            for axis in 0..3 {
                columns.push(SymbolicAffineColumn::Field {
                    dimension: equality.dimension(),
                    support,
                    component: SymbolicFieldComponent::Gradient(axis),
                });
            }
        }
    }
    columns.extend(equality.latent_coefficients().iter().map(|term| {
        SymbolicAffineColumn::SemanticLatent {
            dimension: equality.dimension(),
            latent: term.latent,
        }
    }));
    columns
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalAffineSolverRow {
    pub(crate) row: CanonicalSolverRow,
    pub(crate) sense: CanonicalInequalitySense,
    pub(crate) violation_loss: Option<CanonicalViolationLoss>,
}

impl CanonicalAffineSolverRow {
    pub(crate) fn upper_form_coefficients(
        &self,
        coordinate_layout: CubicFieldCoordinateLayout,
        variable_layout: CubicSolverVariableLayout,
    ) -> Vec<f64> {
        let multiplier = self.sense.upper_form_multiplier();
        let mut coefficients = self.row.coefficients(coordinate_layout, variable_layout);
        coefficients
            .iter_mut()
            .for_each(|coefficient| *coefficient *= multiplier);
        coefficients
    }

    pub(crate) fn upper_form_bound(&self) -> f64 {
        self.sense.upper_form_multiplier() * self.row.target
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftSolverRow {
    pub(crate) row: CanonicalSolverRow,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftObjectiveBlock {
    pub(crate) objective_index: usize,
    pub(crate) canonical_indices: Vec<usize>,
    pub(crate) residuals: Vec<ResidualId>,
    pub(crate) loss: CanonicalSoftLoss,
    pub(crate) precision: Vec<f64>,
    pub(crate) whitening: Vec<f64>,
    pub(crate) inverse_whitening: Vec<f64>,
    pub(crate) covariance_group: Option<GroupId>,
    pub(crate) block_kind: CanonicalSoftResidualBlockKind,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftRecoveryRelation {
    pub(crate) canonical_index: usize,
    pub(crate) provenance: UsageProvenance,
    pub(crate) residual: ResidualId,
    pub(crate) objective_index: usize,
    pub(crate) objective_component: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalSoftRecoveryGraph {
    /// Soft rows are deliberately retained one-for-one until an objective-
    /// preserving compression has its own exact proof and recovery map.
    pub(crate) retained_rows: Vec<usize>,
    pub(crate) relations: Vec<CanonicalSoftRecoveryRelation>,
    rows: Vec<CanonicalSoftSolverRow>,
    objectives: Vec<CanonicalSoftObjectiveBlock>,
}

impl CanonicalSoftRecoveryGraph {
    fn build(rows: &[CanonicalSoftSolverRow], objectives: &[CanonicalSoftObjectiveBlock]) -> Self {
        let relations = objectives
            .iter()
            .flat_map(|objective| {
                objective.canonical_indices.iter().copied().enumerate().map(
                    move |(objective_component, canonical_index)| {
                        let row = &rows[canonical_index].row;
                        CanonicalSoftRecoveryRelation {
                            canonical_index,
                            provenance: row.provenance.clone(),
                            residual: row.residual.clone(),
                            objective_index: objective.objective_index,
                            objective_component,
                        }
                    },
                )
            })
            .collect();
        Self {
            retained_rows: (0..rows.len()).collect(),
            relations,
            rows: rows.to_vec(),
            objectives: objectives.to_vec(),
        }
    }

    fn verifies(
        &self,
        rows: &[CanonicalSoftSolverRow],
        objectives: &[CanonicalSoftObjectiveBlock],
    ) -> bool {
        self.retained_rows.iter().copied().eq(0..rows.len())
            && self.rows == rows
            && self.objectives == objectives
            && self.relations.len() == rows.len()
            && self.relations.iter().enumerate().all(|(index, relation)| {
                let Some(row) = rows.get(relation.canonical_index) else {
                    return false;
                };
                let Some(objective) = objectives.get(relation.objective_index) else {
                    return false;
                };
                relation.canonical_index == index
                    && relation.provenance == row.row.provenance
                    && relation.residual == row.row.residual
                    && objective.objective_index == relation.objective_index
                    && objective
                        .canonical_indices
                        .get(relation.objective_component)
                        == Some(&relation.canonical_index)
                    && objective.residuals.get(relation.objective_component)
                        == Some(&relation.residual)
            })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalCubicSolverForm {
    pub(crate) standard_field_variables: usize,
    pub(crate) quotient_field_variables: usize,
    pub(crate) polynomial_variables: usize,
    pub(crate) semantic_latents: usize,
    pub(crate) standard_field_energy: Vec<f64>,
    pub(crate) quotient_field_energy: Vec<f64>,
    pub(crate) standard_side_conditions: Vec<f64>,
    pub(crate) characteristic_length: f64,
    pub(crate) representation_evidence: CpdEvidence,
    pub(crate) hard_rows: Vec<CanonicalHardSolverRow>,
    pub(crate) hard_recovery: CanonicalHardRecoveryGraph,
    symbolic_hard_rows: Vec<SymbolicAffineRow>,
    pub(crate) affine_rows: Vec<CanonicalAffineSolverRow>,
    pub(crate) soft_rows: Vec<CanonicalSoftSolverRow>,
    pub(crate) soft_objectives: Vec<CanonicalSoftObjectiveBlock>,
    pub(crate) soft_recovery: CanonicalSoftRecoveryGraph,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CanonicalCubicFieldForm {
    standard_field_variables: usize,
    quotient_field_variables: usize,
    polynomial_variables: usize,
    standard_field_energy: Vec<f64>,
    quotient_field_energy: Vec<f64>,
    standard_side_conditions: Vec<f64>,
    characteristic_length: f64,
    representation_evidence: CpdEvidence,
}

impl CanonicalCubicFieldForm {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        standard_field_variables: usize,
        quotient_field_variables: usize,
        polynomial_variables: usize,
        standard_field_energy: Vec<f64>,
        quotient_field_energy: Vec<f64>,
        standard_side_conditions: Vec<f64>,
        characteristic_length: f64,
        representation_evidence: CpdEvidence,
    ) -> Self {
        Self {
            standard_field_variables,
            quotient_field_variables,
            polynomial_variables,
            standard_field_energy,
            quotient_field_energy,
            standard_side_conditions,
            characteristic_length,
            representation_evidence,
        }
    }
}

impl CanonicalCubicSolverForm {
    pub(crate) fn assemble(
        representation: &CubicRepresentation,
        field_form: CanonicalCubicFieldForm,
        problem: &CubicCanonicalProblem,
    ) -> Result<Self, RepresentationFailure> {
        let standard_field_variables = field_form.standard_field_variables;
        let quotient_field_variables = field_form.quotient_field_variables;

        let hard_rows = problem
            .equalities
            .iter()
            .enumerate()
            .map(|(canonical_index, equality)| {
                let provenance = equality.provenance().clone();
                Ok(CanonicalHardSolverRow {
                    row: CanonicalSolverRow::from_canonical(
                        representation,
                        CanonicalSolverRowInput {
                            canonical_index,
                            functional: equality.field().map(|field| field.functional()),
                            latent_coefficients: equality.latent_coefficients(),
                            source_provenances: equality
                                .source_recoveries()
                                .iter()
                                .map(|source| source.provenance.clone())
                                .collect(),
                            provenance,
                            dimension: equality.dimension(),
                            target: equality.target(),
                        },
                    )?,
                    participation: equality.participation(),
                })
            })
            .collect::<Result<Vec<_>, RepresentationFailure>>()?;
        let symbolic_hard_rows = symbolic_hard_rows(problem);
        let hard_recovery = build_hard_recovery_graph(&symbolic_hard_rows);
        debug_assert!(hard_recovery.verifies(&symbolic_hard_rows));

        let affine_rows = problem
            .affine_inequalities
            .iter()
            .enumerate()
            .map(|(canonical_index, inequality)| {
                let provenance = inequality.provenance().clone();
                Ok(CanonicalAffineSolverRow {
                    row: CanonicalSolverRow::from_canonical(
                        representation,
                        CanonicalSolverRowInput {
                            canonical_index,
                            functional: inequality.field().map(|field| field.functional()),
                            latent_coefficients: inequality.latent_coefficients(),
                            source_provenances: inequality.source_provenances().to_vec(),
                            provenance,
                            dimension: inequality.dimension(),
                            target: inequality.bound(),
                        },
                    )?,
                    sense: inequality.sense(),
                    violation_loss: inequality.violation_channel().map(|channel| channel.loss()),
                })
            })
            .collect::<Result<Vec<_>, RepresentationFailure>>()?;

        let soft_rows = problem
            .soft_equalities
            .iter()
            .enumerate()
            .map(|(canonical_index, equality)| {
                let provenance = equality.provenance().clone();
                Ok(CanonicalSoftSolverRow {
                    row: CanonicalSolverRow::from_canonical(
                        representation,
                        CanonicalSolverRowInput {
                            canonical_index,
                            functional: Some(equality.field().functional()),
                            latent_coefficients: &[],
                            source_provenances: vec![provenance.clone()],
                            provenance,
                            dimension: equality.dimension(),
                            target: equality.target(),
                        },
                    )?,
                })
            })
            .collect::<Result<Vec<_>, RepresentationFailure>>()?;

        let mut canonical_offset = 0;
        let soft_objectives: Vec<CanonicalSoftObjectiveBlock> = problem
            .soft_objectives
            .iter()
            .enumerate()
            .map(|(objective_index, objective)| {
                let dimension = objective.residuals().len();
                let canonical_indices =
                    (canonical_offset..canonical_offset + dimension).collect::<Vec<_>>();
                canonical_offset += dimension;
                CanonicalSoftObjectiveBlock {
                    objective_index,
                    canonical_indices,
                    residuals: objective.residuals().to_vec(),
                    loss: objective.loss().clone(),
                    precision: objective.loss().precision_matrix(dimension),
                    whitening: objective.loss().whitening_matrix(dimension),
                    inverse_whitening: objective.loss().inverse_whitening_matrix(dimension),
                    covariance_group: objective.covariance_group().cloned(),
                    block_kind: objective.block_kind().clone(),
                }
            })
            .collect();
        let soft_recovery = CanonicalSoftRecoveryGraph::build(&soft_rows, &soft_objectives);
        debug_assert!(soft_recovery.verifies(&soft_rows, &soft_objectives));

        Ok(Self {
            standard_field_variables,
            quotient_field_variables,
            polynomial_variables: field_form.polynomial_variables,
            semantic_latents: problem.semantic_latents.len(),
            standard_field_energy: field_form.standard_field_energy,
            quotient_field_energy: field_form.quotient_field_energy,
            standard_side_conditions: field_form.standard_side_conditions,
            characteristic_length: field_form.characteristic_length,
            representation_evidence: field_form.representation_evidence,
            hard_rows,
            hard_recovery,
            symbolic_hard_rows,
            affine_rows,
            soft_rows,
            soft_objectives,
            soft_recovery,
        })
    }

    pub(crate) fn variable_layout(
        &self,
        coordinate_layout: CubicFieldCoordinateLayout,
        auxiliary: usize,
    ) -> CubicSolverVariableLayout {
        CubicSolverVariableLayout {
            field: match coordinate_layout {
                CubicFieldCoordinateLayout::Standard => self.standard_field_variables,
                CubicFieldCoordinateLayout::Quotient => self.quotient_field_variables,
            },
            polynomial: self.polynomial_variables,
            semantic_latents: self.semantic_latents,
            auxiliary,
        }
    }

    pub(crate) fn solver_hard_rows(&self) -> impl Iterator<Item = &CanonicalHardSolverRow> {
        self.hard_recovery
            .retained_rows
            .iter()
            .map(|canonical_index| &self.hard_rows[*canonical_index])
    }

    pub(crate) fn verifies_hard_recovery(&self) -> bool {
        self.hard_recovery.verifies(&self.symbolic_hard_rows)
    }

    pub(crate) fn verifies_hard_conflict_witness(&self) -> bool {
        self.hard_recovery
            .verifies_conflict_witness(&self.symbolic_hard_rows)
    }

    pub(crate) fn verifies_soft_recovery(&self) -> bool {
        self.soft_recovery
            .verifies(&self.soft_rows, &self.soft_objectives)
    }

    pub(crate) fn verify_all_source_recovery(
        &self,
        recovered_hard_relations: &[RecoveredHardEquality],
        hard_relation_tolerances: &[CanonicalRelationToleranceEvidence],
        recovered_affine_relations: &[RecoveredAffineRelation],
        recovered_soft_relations: &[RecoveredSoftEquality],
        recovered_soft_objectives: &[RecoveredSoftObjective],
    ) -> AllSourceRecoveryLedger {
        let canonical_hard_relations = self.hard_rows.len()
            + self
                .affine_rows
                .iter()
                .filter(|row| row.violation_loss.is_none())
                .count();
        let canonical_soft_relations = self.soft_rows.len()
            + self
                .affine_rows
                .iter()
                .filter(|row| row.violation_loss.is_some())
                .count();
        let mut participating_sources = self
            .hard_rows
            .iter()
            .flat_map(|row| &row.row.source_provenances)
            .chain(
                self.affine_rows
                    .iter()
                    .flat_map(|row| &row.row.source_provenances),
            )
            .chain(
                self.soft_rows
                    .iter()
                    .flat_map(|row| &row.row.source_provenances),
            )
            .map(|provenance| provenance.source().clone())
            .collect::<Vec<_>>();
        participating_sources.sort();
        participating_sources.dedup();

        let mut recovered_sources = self
            .hard_recovery
            .relations
            .iter()
            .map(|relation| relation.provenance.source().clone())
            .chain(
                recovered_affine_relations
                    .iter()
                    .flat_map(|relation| &relation.provenances)
                    .map(|provenance| provenance.source().clone()),
            )
            .chain(
                recovered_soft_relations
                    .iter()
                    .map(|relation| relation.provenance.source().clone()),
            )
            .collect::<Vec<_>>();
        recovered_sources.sort();
        recovered_sources.dedup();

        let objective_associations_verified = recovered_soft_objectives.len()
            == self.soft_objectives.len()
            && recovered_soft_objectives
                .iter()
                .zip(&self.soft_objectives)
                .all(|(recovered, canonical)| {
                    recovered.canonical_indices == canonical.canonical_indices
                        && recovered.loss == canonical.loss
                        && recovered.covariance_group == canonical.covariance_group
                        && recovered.block_kind == canonical.block_kind
                });
        let every_source_hard_relation_accepted =
            self.hard_recovery.relations.iter().all(|relation| {
                recovered_hard_relations
                    .get(relation.canonical_index)
                    .zip(hard_relation_tolerances.get(relation.canonical_index))
                    .is_some_and(|(recovered, tolerance)| {
                        // Source merging normalizes orientation only. Reject any
                        // future non-unit map until it carries its own scale
                        // evidence; for ±1 the two envelope scales are invariant.
                        let source_tolerance = crate::numerical::EQUALITY_KKT_POLICY_V2
                            .canonical_characteristic_tolerance_multiplier
                            * tolerance.characteristic_scale
                            + crate::numerical::EQUALITY_KKT_POLICY_V2
                                .canonical_relation_reference_tolerance_multiplier
                                * tolerance.relation_reference_scale;
                        relation.relation_to_canonical.abs() == 1.0
                            && (relation.relation_to_canonical * recovered.value - relation.target)
                                .abs()
                                <= source_tolerance
                    })
            });
        let affine_associations_verified = recovered_affine_relations.len()
            == self.affine_rows.len()
            && recovered_affine_relations
                .iter()
                .enumerate()
                .zip(&self.affine_rows)
                .all(|((canonical_index, recovered), canonical)| {
                    recovered.canonical_index == canonical_index
                        && recovered.provenances == canonical.row.source_provenances
                        && recovered.violation_loss == canonical.violation_loss
                        && (recovered.violation_loss.is_some()
                            || recovered.violation <= recovered.tolerance)
                });
        let soft_relations_verified = recovered_soft_relations.len() == self.soft_rows.len()
            && recovered_soft_relations
                .iter()
                .zip(&self.soft_recovery.relations)
                .all(|(recovered, edge)| {
                    recovered.provenance == edge.provenance
                        && edge.canonical_index < self.soft_rows.len()
                });
        let recovery_edges = self.hard_recovery.relations.len()
            + recovered_affine_relations
                .iter()
                .map(|relation| relation.provenances.len())
                .sum::<usize>()
            + self.soft_recovery.relations.len();
        let verified = self.verifies_hard_recovery()
            && self.verifies_soft_recovery()
            && recovered_hard_relations.len() == self.hard_rows.len()
            && hard_relation_tolerances.len() == self.hard_rows.len()
            && every_source_hard_relation_accepted
            && affine_associations_verified
            && soft_relations_verified
            && objective_associations_verified
            && participating_sources == recovered_sources;

        AllSourceRecoveryLedger {
            canonical_hard_relations,
            canonical_soft_relations,
            participating_sources,
            recovered_sources,
            representers: self.representation_evidence.fitting_functional_count,
            solver_relation_rows: self.hard_recovery.retained_rows.len()
                + self.affine_rows.len()
                + self.soft_recovery.retained_rows.len(),
            recovery_edges,
            verified,
        }
    }

    pub(crate) fn field_energy(&self, coordinate_layout: CubicFieldCoordinateLayout) -> &[f64] {
        match coordinate_layout {
            CubicFieldCoordinateLayout::Standard => &self.standard_field_energy,
            CubicFieldCoordinateLayout::Quotient => &self.quotient_field_energy,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{
        CanonicalSoftObjectiveBlock, CanonicalSoftRecoveryGraph, CanonicalSoftSolverRow,
        CanonicalSolverRow, SymbolicAffineColumn, SymbolicAffineRow, SymbolicFieldComponent,
        build_hard_recovery_graph,
    };
    use crate::cubic_equality::{
        CanonicalHardSourceRecovery, CanonicalSoftLoss, CanonicalSoftResidualBlockKind,
        CanonicalSoftResidualMemberKind,
    };
    use crate::functional::{
        DerivedBlockId, DerivedColumnId, DerivedRowId, FunctionalDimension, RelationId, ResidualId,
        SemanticRolePath, SourceId, UsageProvenance,
    };

    fn provenance(source: &str, role: &str) -> UsageProvenance {
        UsageProvenance::new(
            SourceId::new(source),
            None,
            RelationId::new(format!("{source}/relation")),
            ResidualId::new(format!("{source}/residual")),
            SemanticRolePath::new(role),
        )
    }

    #[test]
    fn hard_recovery_requires_exact_row_and_target_reconstruction() {
        let rows = vec![
            SymbolicAffineRow::new(vec![1.0, 0.0], 2.0),
            SymbolicAffineRow::new(vec![0.0, 1.0], 3.0),
            SymbolicAffineRow::new(vec![1.0, 1.0], 5.0),
            SymbolicAffineRow::new(vec![1.0, 1.0], 5.0 + 4.0 * f64::EPSILON),
        ];

        let graph = build_hard_recovery_graph(&rows);

        assert_eq!(graph.retained_rows, vec![0, 1, 3]);
        assert_eq!(graph.rows[0].coefficients, vec![(0, 1.0)]);
        assert_eq!(graph.rows[1].coefficients, vec![(1, 1.0)]);
        assert_eq!(graph.rows[2].coefficients, vec![(0, 1.0), (1, 1.0)]);
        assert_eq!(graph.rows[3].coefficients, vec![(2, 1.0)]);
        assert!(graph.verifies(&rows));
    }

    #[test]
    fn rounded_f64_sum_does_not_count_as_exact_target_reconstruction() {
        let rows = vec![
            SymbolicAffineRow::new(vec![1.0, 0.0], 1.0),
            SymbolicAffineRow::new(vec![0.0, 1.0], 2.0_f64.powi(-54)),
            SymbolicAffineRow::new(vec![1.0, 1.0], 1.0),
        ];

        let graph = build_hard_recovery_graph(&rows);

        assert_eq!(graph.retained_rows, vec![0, 1, 2]);
        assert!(graph.verifies(&rows));
    }

    #[test]
    fn distinct_support_atoms_are_never_approximately_compressed() {
        let rows = vec![
            SymbolicAffineRow::new(vec![1.0, 0.0], 0.0),
            SymbolicAffineRow::new(vec![1.0, f64::MIN_POSITIVE], 0.0),
        ];

        let graph = build_hard_recovery_graph(&rows);

        assert_eq!(graph.retained_rows, vec![0, 1]);
        assert!(graph.verifies(&rows));
    }

    #[test]
    fn symbolic_atoms_keep_physical_dimensions_distinct() {
        let support = [0.0_f64.to_bits(); 3];

        assert_ne!(
            SymbolicAffineColumn::Field {
                dimension: FunctionalDimension::FieldValue,
                support,
                component: SymbolicFieldComponent::Value,
            },
            SymbolicAffineColumn::Field {
                dimension: FunctionalDimension::FieldValuePerLength,
                support,
                component: SymbolicFieldComponent::Value,
            }
        );
    }

    #[test]
    fn recovery_relations_retain_every_source_residual_and_semantic_role() {
        let first = provenance("duplicate-a", "field-value-observation/value");
        let second = provenance("duplicate-b", "additive-field-gauge/point");
        let rows = vec![SymbolicAffineRow {
            coefficients: BTreeMap::from([(0, 1.0)]),
            target: 2.0,
            source_recoveries: vec![
                CanonicalHardSourceRecovery {
                    provenance: first,
                    coefficient: 1.0,
                    target: 2.0,
                },
                CanonicalHardSourceRecovery {
                    provenance: second.clone(),
                    coefficient: -1.0,
                    target: -2.0,
                },
            ],
            solver_constraint: true,
        }];

        let graph = build_hard_recovery_graph(&rows);

        assert_eq!(graph.retained_rows, vec![0]);
        assert_eq!(graph.rows[0].coefficients, vec![(0, 1.0)]);
        assert_eq!(graph.relations[1].coefficients, vec![(0, -1.0)]);
        assert_eq!(graph.relations[1].provenance, second);
        assert_eq!(graph.relations[1].target, -2.0);
        assert_eq!(graph.relations[1].relation_to_canonical, -1.0);
        assert!(graph.verifies(&rows));
    }

    #[test]
    fn hard_conflict_witness_recovers_aliases_and_recomputes_its_margin() {
        let first = provenance("duplicate-a", "field-value-observation/value");
        let alias = provenance("duplicate-b", "additive-field-gauge/point");
        let incompatible = provenance("incompatible", "field-value-observation/value");
        let rows = vec![
            SymbolicAffineRow {
                coefficients: BTreeMap::from([(0, 1.0)]),
                target: 1.0,
                source_recoveries: vec![
                    CanonicalHardSourceRecovery {
                        provenance: first,
                        coefficient: 1.0,
                        target: 1.0,
                    },
                    CanonicalHardSourceRecovery {
                        provenance: alias,
                        coefficient: -1.0,
                        target: -1.0,
                    },
                ],
                solver_constraint: true,
            },
            SymbolicAffineRow {
                coefficients: BTreeMap::from([(0, 1.0)]),
                target: 2.0,
                source_recoveries: vec![CanonicalHardSourceRecovery {
                    provenance: incompatible,
                    coefficient: 1.0,
                    target: 2.0,
                }],
                solver_constraint: true,
            },
        ];

        let graph = build_hard_recovery_graph(&rows);
        let witness = graph.conflict_witness.as_ref().unwrap();
        assert_eq!(
            witness
                .sources
                .iter()
                .map(|source| source.source().as_str())
                .collect::<Vec<_>>(),
            ["duplicate-a", "duplicate-b", "incompatible"]
        );
        assert!(graph.verifies_conflict_witness(&rows));

        let mut damaged = graph.clone();
        damaged.conflict_witness.as_mut().unwrap().separation_margin *= 2.0;
        assert!(!damaged.verifies_conflict_witness(&rows));
    }

    #[test]
    fn exact_dependency_margin_survives_floating_target_cancellation() {
        let mut rows = vec![
            SymbolicAffineRow::new(vec![1.0, 0.0], 1.0e16),
            SymbolicAffineRow::new(vec![0.0, 1.0], 1.0),
            SymbolicAffineRow::new(vec![1.0, 1.0], 1.0e16),
        ];
        for (index, row) in rows.iter_mut().enumerate() {
            row.source_recoveries = vec![CanonicalHardSourceRecovery {
                provenance: provenance(&format!("source-{index}"), "hard/value"),
                coefficient: 1.0,
                target: row.target,
            }];
        }

        let graph = build_hard_recovery_graph(&rows);

        assert_eq!(
            graph.conflict_witness.as_ref().unwrap().separation_margin,
            1.0
        );
        assert!(graph.verifies_conflict_witness(&rows));
    }

    #[test]
    fn soft_recovery_rejects_objective_or_source_association_changes() {
        let soft_provenance = provenance("soft-a", "field-value-observation/value");
        let residual = soft_provenance.residual().clone();
        let rows = vec![CanonicalSoftSolverRow {
            row: CanonicalSolverRow {
                canonical_index: 0,
                response: None,
                latent_coefficients: Vec::new(),
                source_provenances: vec![soft_provenance.clone()],
                derived_block: DerivedBlockId::from_residual(&residual),
                residual: residual.clone(),
                derived_row: DerivedRowId::from_residual(&residual),
                derived_column: Some(DerivedColumnId::from_residual(&residual)),
                dimension: FunctionalDimension::FieldValue,
                target: 2.0,
                provenance: soft_provenance,
            },
        }];
        let objectives = vec![CanonicalSoftObjectiveBlock {
            objective_index: 0,
            canonical_indices: vec![0],
            residuals: vec![residual],
            loss: CanonicalSoftLoss::QuadraticPenalty { weight: 2.0 },
            precision: vec![2.0],
            whitening: vec![2.0_f64.sqrt()],
            inverse_whitening: vec![1.0 / 2.0_f64.sqrt()],
            covariance_group: None,
            block_kind: CanonicalSoftResidualBlockKind::Independent(
                CanonicalSoftResidualMemberKind::FieldValue,
            ),
        }];
        let graph = CanonicalSoftRecoveryGraph::build(&rows, &objectives);
        assert!(graph.verifies(&rows, &objectives));

        let mut changed_objective = objectives.clone();
        changed_objective[0].loss = CanonicalSoftLoss::QuadraticPenalty { weight: 3.0 };
        assert!(!graph.verifies(&rows, &changed_objective));

        let mut changed_source = rows.clone();
        changed_source[0].row.provenance = provenance("soft-b", "field-value-observation/value");
        assert!(!graph.verifies(&changed_source, &objectives));

        let mut changed_target = rows.clone();
        changed_target[0].row.target = 3.0;
        assert!(!graph.verifies(&changed_target, &objectives));

        let mut changed_dimension = rows.clone();
        changed_dimension[0].row.dimension = FunctionalDimension::FieldValuePerLength;
        assert!(!graph.verifies(&changed_dimension, &objectives));
    }
}
