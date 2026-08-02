use std::mem::size_of;

use crate::faer_backend;

pub(crate) const CAPACITY_LIMIT_BYTES: u64 = 8 * 1024 * 1024 * 1024;
pub(crate) const REPORT_FIXED_BYTES: u64 = 512;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CapacityComponent {
    KktDimension,
    Canonical,
    EqualityDense,
    Kkt,
    FactorStorage,
    FaerFactorWorkspace,
    FaerSolveWorkspace,
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
        components: CapacityComponents,
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
pub(crate) struct FaerWorkspaceEvidence {
    pub(crate) factor_bytes: u64,
    pub(crate) factor_alignment: usize,
    pub(crate) solve_bytes: u64,
    pub(crate) solve_alignment: usize,
    pub(crate) peak_bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CapacityComponents {
    pub(crate) canonical_bytes: u64,
    pub(crate) equality_dense_bytes: u64,
    pub(crate) kkt_bytes: u64,
    pub(crate) factor_storage_bytes: u64,
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

pub(crate) fn plan_equality_capacity(
    primal_variables: usize,
    equality_constraints: usize,
) -> Result<EqualityCapacityPlan, CapacityExceededEvidence> {
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
    let equalities = usize_as_u64(equality_constraints, CapacityComponent::Canonical)?;
    let dimension = usize_as_u64(kkt_dimension, CapacityComponent::KktDimension)?;
    let scalar_bytes = size_of::<f64>() as u64;
    let index_bytes = size_of::<usize>() as u64;

    let canonical_scalars = checked_add(CapacityComponent::Canonical, variables, equalities)?;
    let canonical_bytes = checked_mul(
        CapacityComponent::Canonical,
        canonical_scalars,
        scalar_bytes,
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

    let faer_workspace = faer_workspace(kkt_dimension)?;
    let recovery_bytes = checked_mul(
        CapacityComponent::Recovery,
        checked_mul(CapacityComponent::Recovery, dimension, 3)?,
        scalar_bytes,
    )?;
    let report_bytes = checked_add(
        CapacityComponent::Report,
        checked_mul(CapacityComponent::Report, dimension, scalar_bytes)?,
        REPORT_FIXED_BYTES,
    )?;
    let components = CapacityComponents {
        canonical_bytes,
        equality_dense_bytes,
        kkt_bytes,
        factor_storage_bytes,
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
                components,
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
    let factor_alignment = factor.align_bytes();
    let solve_alignment = solve.align_bytes();
    if factor_alignment == 0 {
        return Err(overflow(
            CapacityComponent::FaerFactorWorkspace,
            ArithmeticOperation::FaerStackLayout,
            usize_as_evidence(dimension),
            0,
        ));
    }
    if solve_alignment == 0 {
        return Err(overflow(
            CapacityComponent::FaerSolveWorkspace,
            ArithmeticOperation::FaerStackLayout,
            usize_as_evidence(dimension),
            1,
        ));
    }
    let factor_bytes = usize_as_u64(factor.size_bytes(), CapacityComponent::FaerFactorWorkspace)?;
    let solve_bytes = usize_as_u64(solve.size_bytes(), CapacityComponent::FaerSolveWorkspace)?;

    Ok(FaerWorkspaceEvidence {
        factor_bytes,
        factor_alignment,
        solve_bytes,
        solve_alignment,
        peak_bytes: factor_bytes.max(solve_bytes),
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
        assert_eq!(plan.faer_workspace.factor_bytes, 64);
        assert_eq!(plan.faer_workspace.factor_alignment, 64);
        assert_eq!(plan.faer_workspace.solve_bytes, 64);
        assert!(plan.peak_bytes <= CAPACITY_LIMIT_BYTES);
    }

    #[test]
    fn first_equality_shape_over_eight_gib_is_rejected_before_allocation() {
        assert!(plan_equality_capacity(0, 23_151).is_ok());
        let evidence = plan_equality_capacity(0, 23_152)
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
        assert_eq!(planned_peak_bytes, 8_589_948_160);
        assert_eq!(planned_peak_bytes - evidence.limit_bytes, 13_568);
        assert_eq!(components.faer_workspace.factor_bytes, 12_039_040);
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
}
