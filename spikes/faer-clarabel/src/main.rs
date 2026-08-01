use std::error::Error;
use std::process::Command;

use georbf_backend_probe::{
    probe_clarabel_primal_infeasible, probe_clarabel_qp, probe_clarabel_socp,
    probe_clarabel_unbounded, probe_faer,
};

fn main() -> Result<(), Box<dyn Error>> {
    let faer = probe_faer()?;
    let qp = probe_clarabel_qp()?;
    let socp = probe_clarabel_socp()?;
    let infeasible = probe_clarabel_primal_infeasible()?;
    let unbounded = probe_clarabel_unbounded()?;

    println!("probe.schema=1");
    println!("toolchain.rustc={}", rustc_version()?);
    println!("toolchain.target={}", target_triple());
    println!("faer.version={}", faer.version);
    println!("faer.features={}", faer.features.join(","));
    println!("faer.threads.requested={}", faer.threads.requested);
    println!("faer.threads.actual={}", faer.threads.actual);
    println!(
        "faer.threads.process_global_state_modified={}",
        faer.threads.process_global_state_modified
    );
    println!("faer.kkt.solution={:?}", faer.kkt.solution);
    println!(
        "faer.kkt.normalized_backward_error={:.17e}",
        faer.kkt.normalized_backward_error
    );
    println!(
        "faer.kkt.inertia={}/{}/{}",
        faer.kkt.inertia.positive, faer.kkt.inertia.negative, faer.kkt.inertia.zero
    );
    println!(
        "faer.kkt.scaling.adapter_applied={}",
        faer.kkt.adapter_scaling_applied
    );
    println!(
        "faer.kkt.failure_reason={}",
        faer.failure_reason
            .map_or_else(|| "none".to_owned(), |reason| format!("{reason:?}"))
    );
    println!(
        "faer.rank.qr/svd={}/{}",
        faer.factorizations.col_pivoted_qr_rank, faer.factorizations.svd_rank
    );
    println!(
        "faer.svd.singular_values={:?}",
        faer.factorizations.singular_values
    );
    println!(
        "faer.cholesky.spd/indefinite={}/{}",
        faer.factorizations.cholesky_spd_succeeded,
        faer.factorizations.cholesky_indefinite_rejected
    );
    println!(
        "faer.capacity.state={:?}",
        faer.factorizations.capacity.state
    );
    println!(
        "faer.capacity.first_oversize_square=dimension:{},bytes:{},outcome:{:?}",
        faer.factorizations.capacity.first_oversize_square_dimension,
        faer.factorizations.capacity.first_oversize_square_bytes,
        faer.factorizations.capacity.representable_oversize
    );
    println!(
        "faer.capacity.arithmetic_overflow={:?}",
        faer.factorizations.capacity.arithmetic_overflow
    );

    print_clarabel("qp", &qp);
    print_clarabel("socp", &socp);
    print_clarabel("primal_infeasible", &infeasible.attempt);
    println!(
        "clarabel.primal_infeasible.certificate={:?}",
        infeasible.certificate
    );
    println!(
        "clarabel.primal_infeasible.residual={:.17e}",
        infeasible.certificate_residual
    );
    println!(
        "clarabel.primal_infeasible.cone_violation={:.17e}",
        infeasible.cone_violation
    );
    println!(
        "clarabel.primal_infeasible.separation_margin={:.17e}",
        infeasible.separation_margin
    );
    print_clarabel("unbounded", &unbounded.attempt);
    println!("clarabel.unbounded.certificate={:?}", unbounded.certificate);
    println!(
        "clarabel.unbounded.residual={:.17e}",
        unbounded.certificate_residual
    );
    println!(
        "clarabel.unbounded.cone_violation={:.17e}",
        unbounded.cone_violation
    );
    println!(
        "clarabel.unbounded.descent_margin={:.17e}",
        unbounded.descent_margin
    );

    Ok(())
}

fn print_clarabel(label: &str, evidence: &georbf_backend_probe::ClarabelEvidence) {
    println!("clarabel.{label}.version={}", evidence.version);
    println!("clarabel.{label}.features={}", evidence.features.join(","));
    println!("clarabel.{label}.class={}", evidence.problem_class);
    println!("clarabel.{label}.cones={}", evidence.cones.join(","));
    println!("clarabel.{label}.linear_solver={}", evidence.linear_solver);
    println!("clarabel.{label}.status={}", evidence.termination);
    println!("clarabel.{label}.primal={:?}", evidence.primal);
    println!("clarabel.{label}.dual={:?}", evidence.dual);
    println!("clarabel.{label}.slack={:?}", evidence.slack);
    println!(
        "clarabel.{label}.residual.primal={:.17e}",
        evidence.primal_residual
    );
    println!(
        "clarabel.{label}.residual.dual={:.17e}",
        evidence.dual_residual
    );
    println!(
        "clarabel.{label}.residual.primal_infeasibility={:.17e}",
        evidence.primal_infeasibility_residual
    );
    println!(
        "clarabel.{label}.residual.dual_infeasibility={:.17e}",
        evidence.dual_infeasibility_residual
    );
    println!("clarabel.{label}.gap.abs={:.17e}", evidence.absolute_gap);
    println!("clarabel.{label}.gap.rel={:.17e}", evidence.relative_gap);
    println!("clarabel.{label}.iterations={}", evidence.iterations);
    println!(
        "clarabel.{label}.settings=max_iter:{},max_threads:{},equilibrate:{},refine:{},static_regularization:{},dynamic_regularization:{},direct:{},tol_feas:{:.1e}",
        evidence.settings.max_iterations,
        evidence.settings.max_threads,
        evidence.settings.equilibration,
        evidence.settings.iterative_refinement,
        evidence.settings.static_regularization,
        evidence.settings.dynamic_regularization,
        evidence.settings.direct_solve_method,
        evidence.settings.feasibility_tolerance,
    );
    println!(
        "clarabel.{label}.threads.requested/actual={}/{}",
        evidence.threads.requested, evidence.threads.actual
    );
    println!(
        "clarabel.{label}.threads.process_global_state_modified={}",
        evidence.threads.process_global_state_modified
    );
    println!(
        "clarabel.{label}.scaling.variable={:?}",
        evidence.scaling.variable
    );
    println!(
        "clarabel.{label}.scaling.inverse_variable={:?}",
        evidence.scaling.inverse_variable
    );
    println!(
        "clarabel.{label}.scaling.constraint={:?}",
        evidence.scaling.constraint
    );
    println!(
        "clarabel.{label}.scaling.inverse_constraint={:?}",
        evidence.scaling.inverse_constraint
    );
    println!(
        "clarabel.{label}.scaling.objective={:.17e}",
        evidence.scaling.objective
    );
    println!(
        "clarabel.{label}.failure_reason={}",
        evidence
            .failure_reason
            .map_or_else(|| "none".to_owned(), |reason| format!("{reason:?}"))
    );
}

fn rustc_version() -> Result<String, Box<dyn Error>> {
    let output = Command::new("rustc").arg("--version").output()?;
    if !output.status.success() {
        return Err("rustc --version failed".into());
    }
    Ok(String::from_utf8(output.stdout)?.trim().to_owned())
}

fn target_triple() -> &'static str {
    #[cfg(all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"))]
    return "x86_64-unknown-linux-gnu";
    #[cfg(all(target_arch = "aarch64", target_os = "linux", target_env = "gnu"))]
    return "aarch64-unknown-linux-gnu";
    #[cfg(all(target_arch = "x86_64", target_os = "macos"))]
    return "x86_64-apple-darwin";
    #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
    return "aarch64-apple-darwin";
    #[cfg(all(target_arch = "x86_64", target_os = "windows", target_env = "msvc"))]
    return "x86_64-pc-windows-msvc";
    #[allow(unreachable_code)]
    "unsupported-target"
}
