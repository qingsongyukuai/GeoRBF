use std::mem::size_of;

use crate::faer_backend;

pub(crate) const CAPACITY_LIMIT_BYTES: u64 = 8 * 1024 * 1024 * 1024;
pub(crate) const REPORT_FIXED_BYTES: u64 = 4 * 1024;
const KKT_REPORT_BYTES_PER_DIMENSION: u64 = 8 * size_of::<f64>() as u64;
// Cover the fixed Vec/Box headers and one-term functional storage in the
// canonical/source lowering records and the final HardRelationAssessment.
// Variable caller identity storage is charged separately for each phase.
const SOURCE_LOWERING_FIXED_BYTES_PER_RELATION: u64 = 1024;
const SOURCE_REPORT_FIXED_BYTES_PER_RELATION: u64 = 1024;
// A source/group identity can be cloned into canonical provenance, usage-edge
// provenance, relation/residual identifiers, conflict proof, and final report.
const SOURCE_IDENTIFIER_BYTES_PER_PHASE_MULTIPLIER: u64 = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CapacityComponent {
    KktDimension,
    Canonical,
    EqualityDense,
    Kkt,
    FactorStorage,
    AnalysisAuxiliary,
    SvdRescueStorage,
    FaerFactorWorkspace,
    FaerSolveWorkspace,
    FaerLltWorkspace,
    FaerQrWorkspace,
    FaerRrqrWorkspace,
    FaerSingularValuesWorkspace,
    FaerInertiaWorkspace,
    FaerSvdRescueWorkspace,
    ConvexQpDimension,
    ConvexQpDense,
    ConvexQpCsc,
    ClarabelFactorWorkspace,
    Recovery,
    Report,
    Peak,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ArithmeticOperation {
    Add,
    Multiply,
    FaerStackLayout,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CapacityExceededReason {
    ArithmeticOverflow {
        component: CapacityComponent,
        operation: ArithmeticOperation,
        left: u64,
        right: u64,
    },
    LimitExceeded {
        planned_peak_bytes: u64,
        components: Box<CapacityComponents>,
    },
    ConvexQpLimitExceeded {
        planned_peak_bytes: u64,
        components: Box<ConvexQpCapacityComponents>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CapacityExceededEvidence {
    pub(crate) limit_bytes: u64,
    pub(crate) large_allocation_attempted: bool,
    pub(crate) backend_invocation_attempted: bool,
    pub(crate) reason: CapacityExceededReason,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct WorkspaceLayoutEvidence {
    pub(crate) bytes: u64,
    pub(crate) alignment: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FaerWorkspaceEvidence {
    pub(crate) factor: WorkspaceLayoutEvidence,
    pub(crate) solve: WorkspaceLayoutEvidence,
    pub(crate) llt: WorkspaceLayoutEvidence,
    pub(crate) qr: WorkspaceLayoutEvidence,
    pub(crate) rrqr: WorkspaceLayoutEvidence,
    pub(crate) singular_values: WorkspaceLayoutEvidence,
    pub(crate) inertia: WorkspaceLayoutEvidence,
    pub(crate) svd_rescue: WorkspaceLayoutEvidence,
    pub(crate) peak_bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CapacityComponents {
    pub(crate) canonical_bytes: u64,
    pub(crate) equality_dense_bytes: u64,
    pub(crate) kkt_bytes: u64,
    pub(crate) factor_storage_bytes: u64,
    pub(crate) analysis_auxiliary_bytes: u64,
    pub(crate) svd_rescue_storage_bytes: u64,
    pub(crate) recovery_bytes: u64,
    pub(crate) report_bytes: u64,
    pub(crate) faer_workspace: FaerWorkspaceEvidence,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EqualityCapacityPlan {
    pub(crate) primal_variables: usize,
    pub(crate) equality_constraints: usize,
    pub(crate) kkt_dimension: usize,
    pub(crate) components: CapacityComponents,
    pub(crate) faer_workspace: FaerWorkspaceEvidence,
    pub(crate) peak_bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ConvexQpCapacityComponents {
    pub(crate) canonical_bytes: u64,
    pub(crate) dense_form_bytes: u64,
    pub(crate) csc_realization_bytes: u64,
    pub(crate) factor_workspace_bytes: u64,
    pub(crate) recovery_bytes: u64,
    pub(crate) report_bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ConvexQpCapacityPlan {
    pub(crate) variables: usize,
    pub(crate) constraints: usize,
    pub(crate) canonical_relations: usize,
    pub(crate) components: ConvexQpCapacityComponents,
    pub(crate) peak_bytes: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct EqualityCapacityShape {
    pub(crate) primal_variables: usize,
    pub(crate) equality_constraints: usize,
    pub(crate) canonical_relations: usize,
    pub(crate) source_lowering: SourceStorageShape,
    pub(crate) source_report: SourceStorageShape,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct SourceStorageShape {
    pub(crate) relations: usize,
    pub(crate) identifier_bytes: usize,
}

pub(crate) fn plan_equality_capacity(
    primal_variables: usize,
    equality_constraints: usize,
) -> Result<EqualityCapacityPlan, CapacityExceededEvidence> {
    plan_equality_capacity_for(EqualityCapacityShape {
        primal_variables,
        equality_constraints,
        canonical_relations: equality_constraints,
        source_lowering: SourceStorageShape::default(),
        source_report: SourceStorageShape::default(),
    })
}

pub(crate) fn plan_convex_qp_capacity(
    variables: usize,
    constraints: usize,
    canonical_relations: usize,
) -> Result<ConvexQpCapacityPlan, CapacityExceededEvidence> {
    let variables_u64 = usize_as_u64(variables, CapacityComponent::ConvexQpDimension)?;
    let constraints_u64 = usize_as_u64(constraints, CapacityComponent::ConvexQpDimension)?;
    let relations_u64 = usize_as_u64(canonical_relations, CapacityComponent::Canonical)?;
    let factor_dimension = variables.checked_add(constraints).ok_or_else(|| {
        overflow(
            CapacityComponent::ConvexQpDimension,
            ArithmeticOperation::Add,
            usize_as_evidence(variables),
            usize_as_evidence(constraints),
        )
    })?;
    let factor_dimension_u64 =
        usize_as_u64(factor_dimension, CapacityComponent::ClarabelFactorWorkspace)?;
    let scalar_bytes = size_of::<f64>() as u64;
    let index_bytes = size_of::<usize>() as u64;

    let canonical_scalars =
        checked_add(CapacityComponent::Canonical, variables_u64, relations_u64)?;
    let canonical_bytes = checked_add(
        CapacityComponent::Canonical,
        checked_mul(
            CapacityComponent::Canonical,
            canonical_scalars,
            scalar_bytes,
        )?,
        checked_mul(
            CapacityComponent::Canonical,
            relations_u64,
            SOURCE_LOWERING_FIXED_BYTES_PER_RELATION,
        )?,
    )?;

    let hessian_scalars = checked_mul(
        CapacityComponent::ConvexQpDense,
        variables_u64,
        variables_u64,
    )?;
    let constraint_scalars = checked_mul(
        CapacityComponent::ConvexQpDense,
        constraints_u64,
        variables_u64,
    )?;
    let dense_vector_scalars = checked_sum(
        CapacityComponent::ConvexQpDense,
        &[
            checked_mul(CapacityComponent::ConvexQpDense, variables_u64, 2)?,
            checked_mul(CapacityComponent::ConvexQpDense, constraints_u64, 3)?,
        ],
    )?;
    let dense_form_bytes = checked_mul(
        CapacityComponent::ConvexQpDense,
        checked_sum(
            CapacityComponent::ConvexQpDense,
            &[hessian_scalars, constraint_scalars, dense_vector_scalars],
        )?,
        scalar_bytes,
    )?;

    let upper_hessian_entries = checked_mul(
        CapacityComponent::ConvexQpCsc,
        variables_u64,
        checked_add(CapacityComponent::ConvexQpCsc, variables_u64, 1)?,
    )? / 2;
    let constraint_entries = checked_mul(
        CapacityComponent::ConvexQpCsc,
        constraints_u64,
        variables_u64,
    )?;
    let csc_entries = checked_add(
        CapacityComponent::ConvexQpCsc,
        upper_hessian_entries,
        constraint_entries,
    )?;
    let csc_entry_bytes = checked_add(CapacityComponent::ConvexQpCsc, scalar_bytes, index_bytes)?;
    let column_pointers = checked_mul(
        CapacityComponent::ConvexQpCsc,
        checked_mul(CapacityComponent::ConvexQpCsc, variables_u64, 2)?,
        index_bytes,
    )?;
    let csc_realization_bytes = checked_add(
        CapacityComponent::ConvexQpCsc,
        checked_mul(CapacityComponent::ConvexQpCsc, csc_entries, csc_entry_bytes)?,
        column_pointers,
    )?;

    // Clarabel's qdldl workspace is private. Charge four dense-equivalent
    // KKT blocks plus index vectors, a conservative upper bound for this
    // deliberately dense product seam.
    let factor_square = checked_mul(
        CapacityComponent::ClarabelFactorWorkspace,
        factor_dimension_u64,
        factor_dimension_u64,
    )?;
    let factor_workspace_bytes = checked_add(
        CapacityComponent::ClarabelFactorWorkspace,
        checked_mul(
            CapacityComponent::ClarabelFactorWorkspace,
            checked_mul(CapacityComponent::ClarabelFactorWorkspace, factor_square, 4)?,
            checked_add(
                CapacityComponent::ClarabelFactorWorkspace,
                scalar_bytes,
                index_bytes,
            )?,
        )?,
        checked_mul(
            CapacityComponent::ClarabelFactorWorkspace,
            factor_dimension_u64,
            checked_mul(CapacityComponent::ClarabelFactorWorkspace, index_bytes, 8)?,
        )?,
    )?;

    let recovery_bytes = checked_mul(
        CapacityComponent::Recovery,
        checked_add(
            CapacityComponent::Recovery,
            checked_mul(CapacityComponent::Recovery, variables_u64, 4)?,
            checked_mul(CapacityComponent::Recovery, constraints_u64, 5)?,
        )?,
        scalar_bytes,
    )?;
    let report_bytes = checked_add(
        CapacityComponent::Report,
        REPORT_FIXED_BYTES,
        checked_mul(
            CapacityComponent::Report,
            relations_u64,
            SOURCE_REPORT_FIXED_BYTES_PER_RELATION,
        )?,
    )?;
    let components = ConvexQpCapacityComponents {
        canonical_bytes,
        dense_form_bytes,
        csc_realization_bytes,
        factor_workspace_bytes,
        recovery_bytes,
        report_bytes,
    };
    let peak_bytes = checked_sum(
        CapacityComponent::Peak,
        &[
            canonical_bytes,
            dense_form_bytes,
            csc_realization_bytes,
            factor_workspace_bytes,
            recovery_bytes,
            report_bytes,
        ],
    )?;
    if peak_bytes > CAPACITY_LIMIT_BYTES {
        return Err(CapacityExceededEvidence {
            limit_bytes: CAPACITY_LIMIT_BYTES,
            large_allocation_attempted: false,
            backend_invocation_attempted: false,
            reason: CapacityExceededReason::ConvexQpLimitExceeded {
                planned_peak_bytes: peak_bytes,
                components: Box::new(components),
            },
        });
    }
    Ok(ConvexQpCapacityPlan {
        variables,
        constraints,
        canonical_relations,
        components,
        peak_bytes,
    })
}

pub(crate) fn plan_equality_capacity_for(
    shape: EqualityCapacityShape,
) -> Result<EqualityCapacityPlan, CapacityExceededEvidence> {
    let EqualityCapacityShape {
        primal_variables,
        equality_constraints,
        canonical_relations,
        source_lowering,
        source_report,
    } = shape;
    let kkt_dimension = primal_variables
        .checked_add(equality_constraints)
        .ok_or_else(|| {
            overflow(
                CapacityComponent::KktDimension,
                ArithmeticOperation::Add,
                usize_as_evidence(primal_variables),
                usize_as_evidence(equality_constraints),
            )
        })?;
    let variables = usize_as_u64(primal_variables, CapacityComponent::Canonical)?;
    let canonical_relations = usize_as_u64(canonical_relations, CapacityComponent::Canonical)?;
    let equalities = usize_as_u64(equality_constraints, CapacityComponent::EqualityDense)?;
    let source_lowering_relations =
        usize_as_u64(source_lowering.relations, CapacityComponent::Canonical)?;
    let source_lowering_identifier_bytes = usize_as_u64(
        source_lowering.identifier_bytes,
        CapacityComponent::Canonical,
    )?;
    let source_report_relations = usize_as_u64(source_report.relations, CapacityComponent::Report)?;
    let source_report_identifier_bytes =
        usize_as_u64(source_report.identifier_bytes, CapacityComponent::Report)?;
    let dimension = usize_as_u64(kkt_dimension, CapacityComponent::KktDimension)?;
    let scalar_bytes = size_of::<f64>() as u64;
    let index_bytes = size_of::<usize>() as u64;

    let canonical_scalars =
        checked_add(CapacityComponent::Canonical, variables, canonical_relations)?;
    let canonical_scalar_bytes = checked_mul(
        CapacityComponent::Canonical,
        canonical_scalars,
        scalar_bytes,
    )?;
    let source_lowering_fixed_bytes = checked_mul(
        CapacityComponent::Canonical,
        source_lowering_relations,
        SOURCE_LOWERING_FIXED_BYTES_PER_RELATION,
    )?;
    let source_lowering_identity_bytes = checked_mul(
        CapacityComponent::Canonical,
        source_lowering_identifier_bytes,
        SOURCE_IDENTIFIER_BYTES_PER_PHASE_MULTIPLIER,
    )?;
    let canonical_bytes = checked_sum(
        CapacityComponent::Canonical,
        &[
            canonical_scalar_bytes,
            source_lowering_fixed_bytes,
            source_lowering_identity_bytes,
        ],
    )?;

    let hessian_scalars = checked_mul(CapacityComponent::EqualityDense, variables, variables)?;
    let jacobian_scalars = checked_mul(CapacityComponent::EqualityDense, equalities, variables)?;
    let dense_scalars = checked_sum(
        CapacityComponent::EqualityDense,
        &[hessian_scalars, jacobian_scalars, variables, equalities],
    )?;
    let equality_dense_bytes = checked_mul(
        CapacityComponent::EqualityDense,
        dense_scalars,
        scalar_bytes,
    )?;

    let kkt_scalars = checked_mul(CapacityComponent::Kkt, dimension, dimension)?;
    let kkt_bytes = checked_mul(CapacityComponent::Kkt, kkt_scalars, scalar_bytes)?;

    let factor_vector_bytes = checked_mul(
        CapacityComponent::FactorStorage,
        dimension,
        checked_add(
            CapacityComponent::FactorStorage,
            scalar_bytes,
            checked_mul(CapacityComponent::FactorStorage, 2, index_bytes)?,
        )?,
    )?;
    let factor_storage_bytes = checked_add(
        CapacityComponent::FactorStorage,
        kkt_bytes,
        factor_vector_bytes,
    )?;

    let rrqr_coefficients = checked_mul(
        CapacityComponent::AnalysisAuxiliary,
        dimension,
        usize_as_u64(
            faer_backend::qr_block_size(),
            CapacityComponent::AnalysisAuxiliary,
        )?,
    )?;
    let analysis_auxiliary_scalars = checked_add(
        CapacityComponent::AnalysisAuxiliary,
        rrqr_coefficients,
        checked_mul(CapacityComponent::AnalysisAuxiliary, dimension, 2)?,
    )?;
    let analysis_auxiliary_indices =
        checked_mul(CapacityComponent::AnalysisAuxiliary, dimension, 2)?;
    let analysis_auxiliary_bytes = checked_add(
        CapacityComponent::AnalysisAuxiliary,
        checked_mul(
            CapacityComponent::AnalysisAuxiliary,
            analysis_auxiliary_scalars,
            scalar_bytes,
        )?,
        checked_mul(
            CapacityComponent::AnalysisAuxiliary,
            analysis_auxiliary_indices,
            index_bytes,
        )?,
    )?;
    let svd_rescue_matrix_scalars = checked_mul(
        CapacityComponent::SvdRescueStorage,
        checked_mul(CapacityComponent::SvdRescueStorage, dimension, dimension)?,
        2,
    )?;
    let svd_rescue_storage_bytes = checked_mul(
        CapacityComponent::SvdRescueStorage,
        checked_add(
            CapacityComponent::SvdRescueStorage,
            svd_rescue_matrix_scalars,
            checked_mul(CapacityComponent::SvdRescueStorage, dimension, 2)?,
        )?,
        scalar_bytes,
    )?;

    let faer_workspace = faer_workspace(kkt_dimension)?;
    let recovery_bytes = checked_mul(
        CapacityComponent::Recovery,
        checked_mul(CapacityComponent::Recovery, dimension, 3)?,
        scalar_bytes,
    )?;
    let kkt_report_bytes = checked_mul(
        CapacityComponent::Report,
        dimension,
        KKT_REPORT_BYTES_PER_DIMENSION,
    )?;
    let source_relation_fixed_bytes = checked_mul(
        CapacityComponent::Report,
        source_report_relations,
        SOURCE_REPORT_FIXED_BYTES_PER_RELATION,
    )?;
    let source_identifier_lifecycle_bytes = checked_mul(
        CapacityComponent::Report,
        source_report_identifier_bytes,
        SOURCE_IDENTIFIER_BYTES_PER_PHASE_MULTIPLIER,
    )?;
    let report_bytes = checked_sum(
        CapacityComponent::Report,
        &[
            kkt_report_bytes,
            source_relation_fixed_bytes,
            source_identifier_lifecycle_bytes,
            REPORT_FIXED_BYTES,
        ],
    )?;
    let components = CapacityComponents {
        canonical_bytes,
        equality_dense_bytes,
        kkt_bytes,
        factor_storage_bytes,
        analysis_auxiliary_bytes,
        svd_rescue_storage_bytes,
        recovery_bytes,
        report_bytes,
        faer_workspace: faer_workspace.clone(),
    };
    let peak_bytes = checked_sum(
        CapacityComponent::Peak,
        &[
            canonical_bytes,
            equality_dense_bytes,
            kkt_bytes,
            factor_storage_bytes,
            analysis_auxiliary_bytes,
            svd_rescue_storage_bytes,
            faer_workspace.peak_bytes,
            recovery_bytes,
            report_bytes,
        ],
    )?;

    if peak_bytes > CAPACITY_LIMIT_BYTES {
        return Err(CapacityExceededEvidence {
            limit_bytes: CAPACITY_LIMIT_BYTES,
            large_allocation_attempted: false,
            backend_invocation_attempted: false,
            reason: CapacityExceededReason::LimitExceeded {
                planned_peak_bytes: peak_bytes,
                components: Box::new(components),
            },
        });
    }

    Ok(EqualityCapacityPlan {
        primal_variables,
        equality_constraints,
        kkt_dimension,
        components,
        faer_workspace,
        peak_bytes,
    })
}

fn faer_workspace(dimension: usize) -> Result<FaerWorkspaceEvidence, CapacityExceededEvidence> {
    let factor = faer_backend::factor_workspace_requirement(dimension);
    let solve = faer_backend::solve_workspace_requirement(dimension);
    let llt = faer_backend::llt_workspace_requirement(dimension);
    let qr = faer_backend::qr_workspace_requirement(dimension);
    let rrqr = faer_backend::rrqr_workspace_requirement(dimension);
    let singular_values = faer_backend::singular_values_workspace_requirement(dimension);
    let inertia = faer_backend::inertia_workspace_requirement(dimension);
    let svd_rescue = faer_backend::svd_rescue_workspace_requirement(dimension);
    let factor = stack_layout(factor, CapacityComponent::FaerFactorWorkspace, dimension)?;
    let solve = stack_layout(solve, CapacityComponent::FaerSolveWorkspace, dimension)?;
    let llt = stack_layout(llt, CapacityComponent::FaerLltWorkspace, dimension)?;
    let qr = stack_layout(qr, CapacityComponent::FaerQrWorkspace, dimension)?;
    let rrqr = stack_layout(rrqr, CapacityComponent::FaerRrqrWorkspace, dimension)?;
    let singular_values = stack_layout(
        singular_values,
        CapacityComponent::FaerSingularValuesWorkspace,
        dimension,
    )?;
    let inertia = stack_layout(inertia, CapacityComponent::FaerInertiaWorkspace, dimension)?;
    let svd_rescue = stack_layout(
        svd_rescue,
        CapacityComponent::FaerSvdRescueWorkspace,
        dimension,
    )?;
    let peak_bytes = factor
        .bytes
        .max(solve.bytes)
        .max(llt.bytes)
        .max(qr.bytes)
        .max(rrqr.bytes)
        .max(singular_values.bytes)
        .max(inertia.bytes)
        .max(svd_rescue.bytes);

    Ok(FaerWorkspaceEvidence {
        factor,
        solve,
        llt,
        qr,
        rrqr,
        singular_values,
        inertia,
        svd_rescue,
        peak_bytes,
    })
}

fn stack_layout(
    requirement: faer::dyn_stack::StackReq,
    component: CapacityComponent,
    dimension: usize,
) -> Result<WorkspaceLayoutEvidence, CapacityExceededEvidence> {
    let alignment = requirement.align_bytes();
    if alignment == 0 {
        return Err(overflow(
            component,
            ArithmeticOperation::FaerStackLayout,
            usize_as_evidence(dimension),
            0,
        ));
    }
    Ok(WorkspaceLayoutEvidence {
        bytes: usize_as_u64(requirement.size_bytes(), component)?,
        alignment,
    })
}

fn checked_sum(
    component: CapacityComponent,
    values: &[u64],
) -> Result<u64, CapacityExceededEvidence> {
    values
        .iter()
        .try_fold(0, |sum, value| checked_add(component, sum, *value))
}

fn checked_add(
    component: CapacityComponent,
    left: u64,
    right: u64,
) -> Result<u64, CapacityExceededEvidence> {
    left.checked_add(right)
        .ok_or_else(|| overflow(component, ArithmeticOperation::Add, left, right))
}

fn checked_mul(
    component: CapacityComponent,
    left: u64,
    right: u64,
) -> Result<u64, CapacityExceededEvidence> {
    left.checked_mul(right)
        .ok_or_else(|| overflow(component, ArithmeticOperation::Multiply, left, right))
}

fn usize_as_u64(
    value: usize,
    component: CapacityComponent,
) -> Result<u64, CapacityExceededEvidence> {
    u64::try_from(value).map_err(|_| {
        overflow(
            component,
            ArithmeticOperation::Add,
            usize_as_evidence(value),
            0,
        )
    })
}

fn usize_as_evidence(value: usize) -> u64 {
    u64::try_from(value).unwrap_or(u64::MAX)
}

fn overflow(
    component: CapacityComponent,
    operation: ArithmeticOperation,
    left: u64,
    right: u64,
) -> CapacityExceededEvidence {
    CapacityExceededEvidence {
        limit_bytes: CAPACITY_LIMIT_BYTES,
        large_allocation_attempted: false,
        backend_invocation_attempted: false,
        reason: CapacityExceededReason::ArithmeticOverflow {
            component,
            operation,
            left,
            right,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equality_plan_within_budget_includes_exact_faer_workspace() {
        let plan = plan_equality_capacity(2, 1).expect("small equality plan should fit");

        assert_eq!(plan.kkt_dimension, 3);
        assert_eq!(plan.faer_workspace.factor.bytes, 64);
        assert_eq!(plan.faer_workspace.factor.alignment, 64);
        assert_eq!(plan.faer_workspace.solve.bytes, 64);
        assert!(plan.peak_bytes <= CAPACITY_LIMIT_BYTES);
    }

    #[test]
    fn first_equality_shape_over_eight_gib_is_rejected_before_allocation() {
        let mut accepted = 0;
        let mut rejected = 23_151;
        while rejected - accepted > 1 {
            let candidate = accepted + (rejected - accepted) / 2;
            if plan_equality_capacity(0, candidate).is_ok() {
                accepted = candidate;
            } else {
                rejected = candidate;
            }
        }
        assert_eq!(accepted, 9_841);
        assert_eq!(rejected, 9_842);
        let evidence = plan_equality_capacity(0, rejected)
            .expect_err("the first over-limit shape must be rejected deterministically");

        assert!(!evidence.large_allocation_attempted);
        assert_eq!(evidence.limit_bytes, 8_589_934_592);
        let CapacityExceededReason::LimitExceeded {
            planned_peak_bytes,
            components,
        } = evidence.reason
        else {
            panic!("representable capacity excess must not be reported as overflow");
        };
        assert_eq!(planned_peak_bytes, 8_590_810_592);
        assert_eq!(planned_peak_bytes - evidence.limit_bytes, 876_000);
        assert_eq!(components.faer_workspace.factor.bytes, 5_120_960);
        assert!(components.svd_rescue_storage_bytes > components.factor_storage_bytes);
        assert!(components.faer_workspace.svd_rescue.bytes > 0);
    }

    #[test]
    fn equality_shape_arithmetic_overflow_returns_structured_evidence() {
        let evidence = plan_equality_capacity(usize::MAX, usize::MAX)
            .expect_err("unrepresentable dimensions must be rejected deterministically");

        assert!(!evidence.large_allocation_attempted);
        assert!(!evidence.backend_invocation_attempted);
        assert_eq!(evidence.limit_bytes, CAPACITY_LIMIT_BYTES);
        assert_eq!(
            evidence.reason,
            CapacityExceededReason::ArithmeticOverflow {
                component: CapacityComponent::KktDimension,
                operation: ArithmeticOperation::Add,
                left: u64::MAX,
                right: u64::MAX,
            }
        );
    }

    #[test]
    fn convex_qp_plan_accounts_for_every_materialization_before_backend_entry() {
        let plan = plan_convex_qp_capacity(10, 11, 12)
            .expect("the manufactured dense QP should fit the product limit");

        assert_eq!(plan.variables, 10);
        assert_eq!(plan.constraints, 11);
        assert_eq!(plan.canonical_relations, 12);
        assert!(plan.components.canonical_bytes > 0);
        assert!(plan.components.dense_form_bytes > 0);
        assert!(plan.components.csc_realization_bytes > 0);
        assert!(plan.components.factor_workspace_bytes > 0);
        assert!(plan.components.recovery_bytes > 0);
        assert!(plan.components.report_bytes > 0);
        assert!(plan.peak_bytes <= CAPACITY_LIMIT_BYTES);
    }

    #[test]
    fn convex_qp_over_limit_is_rejected_without_allocation_or_backend_entry() {
        let evidence = plan_convex_qp_capacity(20_000, 20_000, 20_000)
            .expect_err("an oversized dense QP must fail its checked capacity preflight");

        assert!(!evidence.large_allocation_attempted);
        assert!(!evidence.backend_invocation_attempted);
        assert_eq!(evidence.limit_bytes, CAPACITY_LIMIT_BYTES);
        match evidence.reason {
            CapacityExceededReason::ConvexQpLimitExceeded {
                planned_peak_bytes,
                components,
            } => {
                assert!(planned_peak_bytes > CAPACITY_LIMIT_BYTES);
                assert!(components.dense_form_bytes > 0);
                assert!(components.csc_realization_bytes > 0);
                assert!(components.factor_workspace_bytes > 0);
            }
            other => panic!("unexpected capacity failure: {other:?}"),
        }
    }

    #[test]
    fn convex_qp_shape_arithmetic_overflow_returns_structured_evidence() {
        let evidence = plan_convex_qp_capacity(usize::MAX, usize::MAX, usize::MAX)
            .expect_err("unrepresentable QP dimensions must fail deterministically");

        assert!(!evidence.large_allocation_attempted);
        assert!(!evidence.backend_invocation_attempted);
        assert!(matches!(
            evidence.reason,
            CapacityExceededReason::ArithmeticOverflow {
                component: CapacityComponent::ConvexQpDimension,
                operation: ArithmeticOperation::Add,
                ..
            }
        ));
    }
}
