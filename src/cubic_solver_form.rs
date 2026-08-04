//! Solver-independent assembly of canonical Cubic objectives and relations.
//!
//! The representation owns numerical coordinates. This module owns canonical
//! hard/soft semantics and records every source-bearing row once; backend
//! adapters only realize one of the coordinate layouts retained by the form.

use crate::cubic_equality::{
    CanonicalEqualityParticipation, CanonicalInequalitySense, CanonicalSoftLoss,
    CanonicalSoftResidualBlockKind, CanonicalViolationLoss, CpdEvidence, CubicCanonicalProblem,
    CubicFunctionalResponse, CubicRepresentation, RepresentationFailure, SemanticLatentCoefficient,
};
use crate::functional::{
    CanonicalFunctional, DerivedBlockId, DerivedColumnId, DerivedRowId, FunctionalDimension,
    GroupId, ResidualId, UsageProvenance,
};

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
    pub(crate) affine_rows: Vec<CanonicalAffineSolverRow>,
    pub(crate) soft_rows: Vec<CanonicalSoftSolverRow>,
    pub(crate) soft_objectives: Vec<CanonicalSoftObjectiveBlock>,
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
                            source_provenances: vec![provenance.clone()],
                            provenance,
                            dimension: equality.dimension(),
                            target: equality.target(),
                        },
                    )?,
                    participation: equality.participation(),
                })
            })
            .collect::<Result<Vec<_>, RepresentationFailure>>()?;

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
        let soft_objectives = problem
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
            affine_rows,
            soft_rows,
            soft_objectives,
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
        self.hard_rows
            .iter()
            .filter(|row| row.participation == CanonicalEqualityParticipation::SolverConstraint)
    }

    pub(crate) fn field_energy(&self, coordinate_layout: CubicFieldCoordinateLayout) -> &[f64] {
        match coordinate_layout {
            CubicFieldCoordinateLayout::Standard => &self.standard_field_energy,
            CubicFieldCoordinateLayout::Quotient => &self.quotient_field_energy,
        }
    }
}
