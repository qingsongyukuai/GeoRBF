use faer::diag::{DiagMut, DiagRef};
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::lblt::{factor, solve};
use faer::{Conj, MatMut, MatRef, Par};

use crate::capacity::{
    CapacityExceededEvidence, EqualityCapacityPlan, FaerWorkspaceEvidence, plan_equality_capacity,
};

const BACKWARD_ERROR_LIMIT: f64 = 1.0e-11;

#[derive(Debug, Clone, Copy)]
pub(crate) struct EqualityKktSystem<'a> {
    pub(crate) primal_variables: usize,
    pub(crate) equality_constraints: usize,
    pub(crate) hessian: &'a [f64],
    pub(crate) equality_jacobian: &'a [f64],
    pub(crate) stationarity_rhs: &'a [f64],
    pub(crate) equality_rhs: &'a [f64],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BackendFingerprint {
    pub(crate) schema_version: u32,
    pub(crate) crate_name: &'static str,
    pub(crate) crate_version: &'static str,
    pub(crate) features: [&'static str; 2],
    pub(crate) algorithm: &'static str,
    pub(crate) factor_workspace_source: &'static str,
    pub(crate) target_arch: &'static str,
    pub(crate) target_os: &'static str,
    pub(crate) requested_threads: usize,
    pub(crate) actual_threads: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct KktSolveEvidence {
    pub(crate) candidate: Vec<f64>,
    pub(crate) equality_multipliers: Vec<f64>,
    pub(crate) normalized_backward_error: f64,
    pub(crate) capacity: EqualityCapacityPlan,
    pub(crate) workspace: FaerWorkspaceEvidence,
    pub(crate) backend: BackendFingerprint,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum KktInputField {
    Hessian,
    EqualityJacobian,
    StationarityRightHandSide,
    EqualityRightHandSide,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum KktFailure {
    Capacity(CapacityExceededEvidence),
    InvalidLength {
        field: KktInputField,
        expected: usize,
        actual: usize,
    },
    NonFiniteInput {
        field: KktInputField,
        index: usize,
    },
    WorkspaceAllocation {
        phase: &'static str,
        bytes: u64,
        alignment: usize,
    },
    BackendContractViolation {
        normalized_backward_error: f64,
        limit: f64,
    },
}

pub(crate) fn solve_equality_kkt(
    system: &EqualityKktSystem<'_>,
) -> Result<KktSolveEvidence, KktFailure> {
    let capacity = plan_equality_capacity(system.primal_variables, system.equality_constraints)
        .map_err(KktFailure::Capacity)?;
    validate_system(system)?;

    let primal_variables = system.primal_variables;
    let equality_constraints = system.equality_constraints;
    let dimension = capacity.kkt_dimension;
    let matrix_elements = dimension * dimension;
    let mut kkt = vec![0.0; matrix_elements];
    for row in 0..primal_variables {
        for column in 0..primal_variables {
            kkt[row + column * dimension] = system.hessian[row * primal_variables + column];
        }
    }
    for equality in 0..equality_constraints {
        for variable in 0..primal_variables {
            let value = system.equality_jacobian[equality * primal_variables + variable];
            kkt[primal_variables + equality + variable * dimension] = value;
            kkt[variable + (primal_variables + equality) * dimension] = value;
        }
    }

    let original_rhs = system
        .stationarity_rhs
        .iter()
        .chain(system.equality_rhs)
        .copied()
        .collect::<Vec<_>>();
    let mut solution = original_rhs.clone();
    let mut factors = kkt.clone();
    let mut subdiagonal = vec![0.0; dimension];
    let mut permutation = vec![0usize; dimension];
    let mut inverse_permutation = vec![0usize; dimension];

    let factor_requirement =
        factor::cholesky_in_place_scratch::<usize, f64>(dimension, Par::Seq, Default::default());
    let mut factor_memory =
        MemBuffer::try_new(factor_requirement).map_err(|_| KktFailure::WorkspaceAllocation {
            phase: "factor",
            bytes: capacity.faer_workspace.factor_bytes,
            alignment: capacity.faer_workspace.factor_alignment,
        })?;
    let (_, permutation) = factor::cholesky_in_place(
        MatMut::from_column_major_slice_mut(&mut factors, dimension, dimension),
        DiagMut::from_slice_mut(&mut subdiagonal),
        &mut permutation,
        &mut inverse_permutation,
        Par::Seq,
        MemStack::new(&mut factor_memory),
        Default::default(),
    );
    drop(factor_memory);

    let solve_requirement = solve::solve_in_place_scratch::<usize, f64>(dimension, 1, Par::Seq);
    let mut solve_memory =
        MemBuffer::try_new(solve_requirement).map_err(|_| KktFailure::WorkspaceAllocation {
            phase: "solve",
            bytes: capacity.faer_workspace.solve_bytes,
            alignment: capacity.faer_workspace.solve_alignment,
        })?;
    solve::solve_in_place_with_conj(
        MatRef::from_column_major_slice(&factors, dimension, dimension),
        MatRef::from_column_major_slice(&factors, dimension, dimension).diagonal(),
        DiagRef::from_slice(&subdiagonal),
        Conj::No,
        permutation,
        MatMut::from_column_major_slice_mut(&mut solution, dimension, 1),
        Par::Seq,
        MemStack::new(&mut solve_memory),
    );

    let normalized_backward_error =
        normalized_backward_error(&kkt, dimension, &solution, &original_rhs);
    if !normalized_backward_error.is_finite()
        || normalized_backward_error > BACKWARD_ERROR_LIMIT
        || solution.iter().any(|value| !value.is_finite())
    {
        return Err(KktFailure::BackendContractViolation {
            normalized_backward_error,
            limit: BACKWARD_ERROR_LIMIT,
        });
    }

    Ok(KktSolveEvidence {
        candidate: solution[..primal_variables].to_vec(),
        equality_multipliers: solution[primal_variables..].to_vec(),
        normalized_backward_error,
        workspace: capacity.faer_workspace.clone(),
        backend: backend_fingerprint(),
        capacity,
    })
}

fn validate_system(system: &EqualityKktSystem<'_>) -> Result<(), KktFailure> {
    let expected_hessian = system
        .primal_variables
        .checked_mul(system.primal_variables)
        .unwrap_or(usize::MAX);
    let expected_jacobian = system
        .equality_constraints
        .checked_mul(system.primal_variables)
        .unwrap_or(usize::MAX);
    validate_slice(KktInputField::Hessian, system.hessian, expected_hessian)?;
    validate_slice(
        KktInputField::EqualityJacobian,
        system.equality_jacobian,
        expected_jacobian,
    )?;
    validate_slice(
        KktInputField::StationarityRightHandSide,
        system.stationarity_rhs,
        system.primal_variables,
    )?;
    validate_slice(
        KktInputField::EqualityRightHandSide,
        system.equality_rhs,
        system.equality_constraints,
    )
}

fn validate_slice(field: KktInputField, values: &[f64], expected: usize) -> Result<(), KktFailure> {
    if values.len() != expected {
        return Err(KktFailure::InvalidLength {
            field,
            expected,
            actual: values.len(),
        });
    }
    if let Some((index, _)) = values
        .iter()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(KktFailure::NonFiniteInput { field, index });
    }
    Ok(())
}

fn normalized_backward_error(
    matrix: &[f64],
    dimension: usize,
    solution: &[f64],
    rhs: &[f64],
) -> f64 {
    let residual_norm = (0..dimension)
        .map(|row| {
            let product = (0..dimension)
                .map(|column| matrix[row + column * dimension] * solution[column])
                .sum::<f64>();
            (product - rhs[row]).abs()
        })
        .fold(0.0_f64, f64::max);
    let matrix_norm = (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix[row + column * dimension].abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let solution_norm = solution
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let rhs_norm = rhs.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    residual_norm / (matrix_norm * solution_norm + rhs_norm)
}

fn backend_fingerprint() -> BackendFingerprint {
    BackendFingerprint {
        schema_version: 1,
        crate_name: "faer",
        crate_version: "0.24.4",
        features: ["linalg", "std"],
        algorithm: "LBLT Bunch-Kaufman",
        factor_workspace_source: "faer::linalg::cholesky::lblt::factor::cholesky_in_place_scratch",
        target_arch: std::env::consts::ARCH,
        target_os: std::env::consts::OS,
        requested_threads: 1,
        actual_threads: Par::Seq.degree(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::capacity::REPORT_FIXED_BYTES;

    #[test]
    fn symmetric_equality_kkt_returns_verified_candidate_and_backend_evidence() {
        let system = EqualityKktSystem {
            primal_variables: 2,
            equality_constraints: 1,
            hessian: &[2.0, 0.0, 0.0, 2.0],
            equality_jacobian: &[1.0, 1.0],
            stationarity_rhs: &[2.0, 0.0],
            equality_rhs: &[0.0],
        };

        let evidence = solve_equality_kkt(&system).expect("manufactured KKT should solve");

        assert_eq!(evidence.candidate, vec![0.5, -0.5]);
        assert_eq!(evidence.equality_multipliers, vec![1.0]);
        assert!(evidence.normalized_backward_error <= 1.0e-11);
        assert_eq!(evidence.capacity.kkt_dimension, 3);
        assert_eq!(evidence.workspace.factor_bytes, 64);
        assert_eq!(evidence.workspace.solve_bytes, 64);
        assert_eq!(evidence.backend.schema_version, 1);
        assert_eq!(evidence.backend.crate_name, "faer");
        assert_eq!(evidence.backend.crate_version, "0.24.4");
        assert_eq!(evidence.backend.features, ["linalg", "std"]);
        assert_eq!(evidence.backend.algorithm, "LBLT Bunch-Kaufman");
        assert_eq!(evidence.backend.requested_threads, 1);
        assert_eq!(evidence.backend.actual_threads, 1);
    }

    #[test]
    fn capacity_report_envelope_covers_the_fixed_kkt_evidence() {
        assert!(std::mem::size_of::<KktSolveEvidence>() as u64 <= REPORT_FIXED_BYTES);
    }
}
